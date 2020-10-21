% =========================================================================
% FORMAT function [Z, R] = nk_AdjForCovarsUsingLocalModels (SY, S, Sid, TY, T, Tid, ...
%                                                tolvec, numsel, wins, act)
% =========================================================================
% performs one-2-multi matching of T to S for normalising TY to SY using
% local models of SY (either using local z-normalisation models or partial
% correlation models
% 
% Inputs:
% -------
% SY        : Source Data Matrix
% S         : Source Covariate Matrix => [n rows x c columns], n = subjects, 
%               c = covariates
% Sid       : IDs of Source Subjects
% TY        : Target Data Matrix
% T         : Target Covariate Matrix => [m rows x c columns], m = subjects
%               c = covariates
% Tid       : IDs of Target Matrix
% tolvec    : Covariate tolerance vecor determing the 'width' of the local 
%               reference population (e.g +/- 5 yrs) 
% numsel    : Max number of subjects to be selected as reference
%               (empty = all possible subjects)
% wins      : Winsorization cutoff for z-normalization
% act       : action string ('znorm', 'partcorr')
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Nikolaos Koutsouleris, 01/2015
function [Z, R] = nk_AdjForCovarsUsingLocalModels(SY, S, Sid, TY, T, Tid, ...
                                                tolvec, numsel, wins, act)

if isempty(TY) || isempty(T) || isempty(Tid)
    modestr = 'loo';
    T = S; Tid = Sid; TY = SY;
else
    modestr = 'ext';
end

[nT, mT] = size(T);
[nS, mS] = size(S);

if ~exist('tolvec','var') || isempty(tolvec)
    tolvec = zeros(1,mT);
    for i=1:mT
        tolvec(i) = nk_Range(T(:,i));
    end
end

nV = size(tolvec,2);

if nV ~= mT && nV ~= mS
    error('Check number of columns across S, T and tolvec')
end

if ~exist('numsel','var'), numsel = []; end
if ~exist('wins','var'), wins = 10; end
Z = zeros(size(TY));
R = cell(nT,1);

ind_nz = any(tolvec,1);

% Loop through all target subjects
fprintf('\nWorking on normalisation')
for i=1:nT
    
    
    % Get current subject
    Ti = T(i,:); 
    
    % Work in LOO mode if needed
    if strcmp(modestr,'loo')
        fprintf('\nLOO mode: ')
        tS = S; tSid = Sid; tSY = SY;
        tS(i,:) = []; tSid(i) = []; tSY(i,:) = [];
    else
        fprintf('\n')
    end
    
    % Build lower and upper tolerance vector
    Ui  = Ti + tolvec;
    Li  = Ti - tolvec;
    
    fprintf('Building ref pop for %s: ',Tid{i});
    fprintf('Upper:'); fprintf(' %3.0f',Ui); fprintf('; ');
    fprintf('Lower:'); fprintf(' %3.0f',Li); fprintf('.');
    
    % Apply tolerance vectors to source matrix
    UiS = bsxfun(@le,tS,Ui);
    LiS = bsxfun(@ge,tS,Li);
    
    % Combine UiS & LiS
    iS  = UiS & LiS;
    
    siS = sum(iS,2); fiS = siS == nV; tiS = sum(fiS);
    fprintf('\t==> %g matches ...', tiS);
    
    R{i}.target_ID = Tid{i};
    R{i}.ref_ID = tSid(fiS);
    R{i}.ref_N = tiS;
    R{i}.ref_rows = find(fiS);
    R{i}.ref_S = tS(fiS,:);
    
    rind = 1:tiS;
    if ~isempty(numsel)
        if tiS > numsel(2);
            % Randomly select numsel(2) subjects from reference population
            rind = randperm(tiS, numsel(2));
            R{i}.ref_ID_sel = tSid{fiS(rind)};
            R{i}.ref_rows_sel = R{i}.ref_rows(rind);
            R{i}.ref_S_sel = tS(R{i}.ref_rows(rind),:);
            
        end
    end
    
    % Get Reference Data and Covariate matrices & normalize data
    Yi = tSY(R{i}.ref_rows(rind),:);
    Si = tS(R{i}.ref_rows(rind),:);
    
    switch act
        
        case 'znorm'
            
            fprintf(' Z-normalising data using median ...')
            
            mYi = median(Yi); 
            sYi = std(Yi);

            Zi = (TY(i,:) - mYi)./sYi;
            Zi( isinf(Zi) | isnan(Zi) ) = 0;

            if ~isempty(wins)
                Zi( Zi > wins ) = wins;
                Zi( Zi < -1*wins ) = -1*wins;
            end
            
        case 'partcorr'
            
            fprintf(' Residualising data (+intercept) using ref pop ...')

            % Compute beta(s) from reference data and covariate matrices 
            [~,beta] = nk_PartialCorrelations(Si(:,ind_nz), Yi);
            
            % Apply beta(s) to given target subject
            Zi = nk_PartialCorrelations(Ti(ind_nz), TY(i,:), beta);
            %R{i}.beta_hist = hist(beta();
    end
    Z(i,:) = Zi;
    fprintf(' Done.')
end
fprintf('\nCompleted!\n')
    
end

% =========================================================================
% FORMAT [Y, beta] = nk_PartialCorrelations(G, Y, beta, nointercept, revertflag)
% =========================================================================
% Remove nuisance effects G from Y (optionally, using a predefined beta)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Nikolaos Koutsouleris, 07 / 2011
function [Y, beta] = nk_PartialCorrelations(G, Y, beta, nointercept, revertflag, subgroup)


if (~exist('nointercept','var') || isempty(nointercept) || ~nointercept ) && ~isequal(G, ones(size(G,1),1))
    %Create intercept vecotr
    interceptflag = true;
    intercept = ones(size(G,1),1);
    %Check if intercept is already included in matrix to avoid double
    %intercept removal
    if exist('beta','var') && ~isempty(beta)
        if size(beta,1) == size(G,2), interceptflag = false; end
    end
else
    interceptflag = false;
end

if interceptflag
    %fprintf(' ... adding intercept to covariate matrix')
    G = [intercept G];
end

if ~exist('beta','var') || isempty(beta), 
    if ~exist('subgroup','var') || isempty(subgroup)
        % Compute beta from entire population
        beta = pinv(G) * Y; 
    else
        % Compute beta from a subgroup of observations
        nG = size(G,2);
        beta = zeros(size(G,2), size(Y,2));
        for i= 1 : nG
            beta(i,:) = pinv(G(subgroup(:,i),i)) * Y(subgroup(:,i),:);
        end
    end
end
if ~exist('revertflag','var') || isempty(revertflag) || ~revertflag
    Y = Y - G * beta;
else
    Y = Y + G * beta;
end

end