% =========================================================================
% function [P, Perf, d, crit] = nk_MultiDecideMaxWins(H, Y, Classes,
%                                                       maxfunc, weightflag)
% =========================================================================
% 
% INPUTS:
% -------
% H :                   Ensemble consisting of binary one-vs-one dichotomizer
%                       prediction for given observations
% Y :                   Multi-group label vector
% Classes :             Dichotomization vector with 1 -> number of dichotomizers in
%                       ensemble
% maxfunc :             Decision function, with 1 = sum, 2 = mean, 3 = product, 
%                       4 = majority vote, 5 = median
% weightflag :          0 = no weighting of decision values by number of
%                       predictions per dichotomization class
% 
% OUTPUTS:
% --------
% P :                   Multi-group predictions
% Perf :                Multi-group prediction performance
% d :                   decision matrix
% crit :                selected first-rank (max) decision score from d
%
% COMMENTS:
% ---------
% This function performs One-Vs-One-Max-Wins multi-group classification
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Nikolaos Koutsouleris, 02/2011

function [P, Perf, d, crit] = nk_MultiDecideMaxWinsNegTiesOpt(H, Y, Classes, maxfunc, weightflag)

global RAND

if iscell(H)
    [iy,jy,nclass]  = size(H);
    d               = cell(iy,jy);
    P               = cell(iy,jy);
    crit            = cell(iy,jy);
    Perf            = zeros(iy,jy);
else
    iy = 1; jy=1;
    nclass          = max(Classes);
end

switch maxfunc
    case 1
        maxfnc= 'sum';
    case 2
        maxfnc= 'meannonzero';
    case 3
        maxfnc= 'prodnonzero';
    case 4
        maxfnc= 'majvote';
    case 5
        maxfnc= 'mediannonzero';
end

% Define processing mode:
% Decompose = 1 => One-vs-One
% Decompose = 2 => One-vs-All
if ~isfield(RAND,'Decompose'), 
    codefxname = 'nk_OneVsOne'; assgnfxname = 'AssignOneVsOne';
else
    switch RAND.Decompose
        case 1
             codefxname = 'nk_OneVsOne'; assgnfxname = 'AssignOneVsOne';
        case 2
             codefxname = 'nk_OneVsAll'; assgnfxname = 'AssignOneVsAll';
    end
   
end

for k=1:iy % Loop through partitions (in case iscell(H) = true)
    
    for l=1:jy % Loop through folds (in case iscell(H) = true)
        
        if iscell(H)
            Hkl = H{k,l};
            Ckl = Classes{k,l};
        else
            Hkl = H;
            Ckl = Classes;
        end
        
        % Construct coding matrix
        M = feval(codefxname,Ckl,nclass);
        nsubj       = size(Hkl,1);
        ngroups     = size(M,1);
        d           = zeros(nsubj,ngroups);
        
        Weights = zeros(1,ngroups);
        if iscell(H), tH= []; for curclass=1:nclass, tH = [tH H{k,l,curclass}]; end; else tH = H; end
        
        td = zeros(nsubj,nclass);
        
        for j=1:nclass
            % Perform either one-vs-one or one-vs-all preprocessing
            [td, Weights] = feval(assgnfxname, td, j, M, Classes, Weights, tH, maxfnc);
        end
        
        % Weight decision structure according to number of dichotomizers in
        % each binary comparison
        if weightflag
            [~,mXI] = max(Weights); Weights = Weights./Weights(mXI); 
            td = td./repmat(Weights,nsubj,1);
        end
        
        % Get negative ties
        xtd = sum(sign(td),2) == -1*size(td,2);
        xmin = zeros(size(xtd));
        %xminl = false(size(td));
        % Get the index to the binary classifier with the maximum value
        [~, xmin(xtd)] = max(td(xtd,:),[],2);
        
        % Create logical array of overruling classifiers
        % Overrule 
        
        
        % Compute multi-group prediction according to maximum score in td
        % score difference between groups might also be taken into account
        % in future versions
        if iscell(H)
            d{k,l} = td; 
            [crit{k,l},P{k,l}] = max(td,[],2);
            if isempty(Y), Y = ones(numel(P{k,l}),1); end
            errs = P{k,l} ~= Y;
        else
            d =td;
            [crit,P] = max(td,[],2);
            if isempty(Y), Y = ones(numel(P),1); end
            errs = P ~= Y;
        end
        
        % Compute multi-group performance as accuracy
        Perf(k,l) = (1-sum(errs)/length(errs))*100;
    end
end

return

% =========================================================================
%%%% ASSIGNMENT FUNCTION : ONE-VS-ONE %%%%
function [td, Weights] = AssignOneVsOne(td, j, M, Classes, Weights, tH, maxfunc)

 % Get binary classifiers for current dichotomizer
indC        = Classes==j;    % Index to current dichotomizers
weightC     = sum(indC);     % Weight of current dichotomizer
tHX         = tH(:,indC);    % Scores of current dichotomizer
mJ          = M(:,indC);     % Extract current dichotomizers from multi-group coding matrix
indG        = find(mJ(:,1)); % Indices to groups
                   
% Produce some decision value according to the following
% functions:
ind1        = tHX > 0;          % Dichotomizers voting for +1 group
ind2        = tHX < 0;          % Dichotomizers voting for -1 group
tHX1        = zeros(size(tHX));
tHX2        = tHX1;
tHX1(ind1)  = abs(tHX(ind1));   % Get respective absolute scores for +1 group
tHX2(ind2)  = abs(tHX(ind2));   % Get respective absolute scores for -1 group
Weights(indG) = Weights(indG) + weightC; % Update weights

if size(tHX,2) == 1, 
    if strcmp(maxfunc,'majvote')
        td(:, indG(1)) = td(:, indG(1)) + ind1; % Majority voting for +1 group
        td(:, indG(2)) = td(:, indG(2)) + ind2; % Majority voting for -1 group
    else
        td(:, indG(1)) = td(:, indG(1)) + tHX1; % Some other function for score evaluation of +1 group
        td(:, indG(2)) = td(:, indG(2)) + tHX2; % Some other function for score evaluation of -1 group
    end
else
    % Update Decision structure td
    if strcmp(maxfunc,'majvote')
        td(:, indG(1)) = td(:, indG(1)) + sum(ind1, 2); % Majority voting for +1 group
        td(:, indG(2)) = td(:, indG(2)) + sum(ind2, 2); % Majority voting for -1 group
    else
        td(:, indG(1)) = td(:, indG(1)) + feval(maxfunc, tHX1, 2); % Some other function for score evaluation of +1 group
        td(:, indG(2)) = td(:, indG(2)) + feval(maxfunc, tHX2, 2); % Some other function for score evaluation of -1 group
    end
end
return

% =========================================================================
%%%% ASSIGNMENT FUNCTION : ONE-VS-ALL %%%%
function [td, Weights] = AssignOneVsAll(td, j, M, Classes, Weights, tH, maxfunc)

% Get binary classifiers for current dichotomizer
indC        = Classes==j;    % Index to current dichotomizers
tHX         = tH(:,indC);    % Scores of current dichotomizer
Weights(j)  = sum(indC);
% ind1        = tHX > 0;
% tHX1        = zeros(size(tHX));
% tHX1(ind1)  = tHX(ind1);
% if isempty(tHX1), return; end

% Update Decision structure td
if strcmp(maxfunc,'majvote')
    td(:, j) = sum(sign(tHX), 2); % Majority voting for +1 group
else                    
    td(:, j) = feval(maxfunc, tHX, 2); % Some other function for score evaluation of +1 group
end

return

% _________________________________________________________________________
function M = meannonzero(X, dim)

M = sum(X, dim)./sum(X~=0, dim);

return
% _________________________________________________________________________
function M = prodnonzero(X, dim)

M = sum(X, dim)./sum(X~=0, dim);

return
% _________________________________________________________________________
function M = mediannonzero(X, dim)

X(X==0) = NaN;
M = nm_nanmedian(X,dim);

return

