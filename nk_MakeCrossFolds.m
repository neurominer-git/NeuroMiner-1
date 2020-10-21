% =========================================================================
% cv = nk_MakeCrossFolds(label, RAND, ...
%                        decomposeflag, modeflag, groups, groupnames, ...
%                        oldcv, appendfl, auto_adjust)
% =========================================================================
% 
% The function produces an outer/inner CV structure using random
% resampling. To assess the correct indices for each inner CV
% training/testing fold the following syntax has to be used:
%
% EXAMPLES:
% cv.TrainInd{1,1}(cv.class{1,1}{1}.TestInd{1,1}))
%
% This will give you the row indices of the inner CV data fold [1,1] 
% in the binary class{1} for the outer cross-validation fold [1,1].
% 
% The outer CV data fold [1,1] of binary classifier 1 is accessed via:
% cv.TestInd{1,1}(cv.classnew{1,1}{1}.ind)
%
% If you want to access one inner fold [k,l] of a binary classifier b of the
% current outer CV index [i,j] you use the expression:
% Y{i,j}.inY{k,l}(dat.class{i,j}{b}.ind,:)
%
% First generate outer CV folds. Whole data is split intop training folds
% (used for model generation / learning and validation fold, that are
% completely unseen by the learning algorithm)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NeuroMiner 1.0, (c) Nikolaos Koutsouleris, 08/2018

function cv = nk_MakeCrossFolds(label, RAND, modeflag, groups, groupnames, oldcv, appendfl, auto_adjust)

OutPerms        = RAND.OuterPerm;
OutFold         = RAND.OuterFold;
InPerms         = RAND.InnerPerm;
InFold          = RAND.InnerFold;
   
% One vs One, One vs All, Multigroup 
decomposeflag   = RAND.Decompose;
CV2LCO = []; if isfield(RAND,'CV2LCO'), CV2LCO = RAND.CV2LCO; end
CV1LCO = []; if isfield(RAND,'CV1LCO'), CV1LCO = RAND.CV1LCO; end
if isfield(RAND,'CV2FRAME') && RAND.CV2FRAME == 4 , cv2frame = true; else, cv2frame = false; end
if ~exist('auto_adjust', 'var'), auto_adjust = false; end

if isempty(CV2LCO)
    if OutFold == -1 || OutFold == numel(label)
        CV2LOO = true;
        OutPerms = 1;
    else
        CV2LOO = false;
    end
end

if InFold == -1 || InFold == numel(label)
    CV1LOO = true;
    InPerms = 1;
else
    CV1LOO = false;
end

% Number of class labels
if ~exist('modeflag','var') || isempty(modeflag), modeflag = 'classification'; end
if ~exist('appendfl','var') || isempty(appendfl), appendfl = false; end

% if ~CV2LOO && ~CV1LOO
%     constraintfl = nk_input('Constrain stratification according to some grouping variable', ...
%         0,'yes|no',[1,0],1);
% else
    constraintfl = false;
%end

if constraintfl
    n_subjects_all = size(label,1);
    Constraint = nk_input('Specify grouping vector',0,'r',[],[n_subjects_all 1]);
else
    Constraint = [];
end
ulb = unique(label,'rows');
if any(~isfinite(ulb)), 
    NaNflag = true; ind = logical(sum(isfinite(ulb),2));
    ulb = ulb(ind,:);
else
    NaNflag = false;
end

if strcmp(modeflag,'classification') && decomposeflag ~= 9
    nclass = numel(ulb);
    if ~exist('oldcv','var') || isempty(oldcv)
        if ~exist('groupnames','var') || isempty(groupnames)
            g = cell(nclass,1);
            groupnameflag = nk_input('Define classifier descriptions',0,'yes|no',[1,0]);

            if groupnameflag
                for i=1:nclass
                    g{i} = nk_input(['Group #' num2str(i) ' name'],0,'s');
                end
            else
                g{nclass} = [];
            end
        else
            g = groupnames;
        end
    else
        b = [];
        nbincomp =  nclass*(nclass-1)/2;
        for i=1:nbincomp
            groupsstr = regexp(oldcv.class{1,1}{i}.groupdesc, ' vs ','split');
            b = [b groupsstr]; 
        end
        [dum, i] = unique(b);
        g = b(sort(i));
    end
    Label = label;
else
    if exist('groups','var')
        if ~isempty(groups), 
            Label = groups; 
        else
            Label = ones(size(label,1),1);
        end;
    else
        Label = ones(length(label),1);
    end
    indnan = logical(sum(isnan(label),2));
    if any(indnan), Label(indnan,:)=NaN; end
end

if ~isempty(CV2LCO)
    fprintf('\nGenerating outer (CV2) independent group / site partitioning.')
elseif CV2LOO
    fprintf('\nGenerating outer (CV2) LOO partitioning.')
else
    fprintf('\nGenerating outer (CV2) crossvalidation partitioning.')
    fprintf('\nK-fold: %g, Perms: %g', OutFold, OutPerms)
end

switch appendfl
    case {0, 1}
        if ~isempty(CV2LCO)
            % Independent group / site validation: Enables e.g. the
            % leave-center-out validation of prediction systems
            cv = nk_INDEPpartition(CV2LCO.ind, Label, cv2frame, OutPerms);
        else
            switch CV2LOO
                case false
                    % This is stratified k-fold cross-validation
                    cv = nk_CVpartition2(OutPerms, OutFold, Label, Constraint);
                    if ~isstruct(cv)
                        OutFold = cv; cv = nk_CVpartition2(OutPerms, OutFold, Label, Constraint);
                        InFold = OutFold - 1;
                    end
                case true
                    % This is LOO cross-validation
                    cv = nk_LOOpartition(Label);
            end
        end
        if isempty(cv), return, end
    case {2,3}
        cv = oldcv;
end
% cv is a struct containing TrainInd/TestInd vectors (row indices) and
% Train/TestFold vectors assigning the subjects to each of the defined fold
% per permutation

% Then generate inner CV folds for each outer training fold 
% This will split the training samples further into training / CV
% folds) used for model generation.

[ix,jx] = size(cv.TrainInd);

if isfield(RAND,'Eq') && RAND.Eq.enabled
    Eq = RAND.Eq;
    Eq.AddRemoved2Test = RAND.Eq.addremoved2test;
    if isfield(RAND.Eq,'maxcount'), 
        Eq.MaxCount = RAND.Eq.maxcount;
    else
        Eq.MaxCount = ceil(numel(Eq.Covar)/40);
    end
    if isfield(RAND.Eq,'mincount'), 
        Eq.MinCount = RAND.Eq.mincount;
    else
        Eq.MinCount = ceil(numel(Eq.Covar)/40);
    end
    if isfield(RAND.Eq,'bincount'), 
        Eq.BinCount = RAND.Eq.bincount;
    else
        Eq.BinCount = 7;
    end
else
    Eq = [];
end

for i=1:ix
    
    for j=1:jx
        
        if ~isempty(Eq)
            Eq.Covar = RAND.Eq.Covar(cv.TrainInd{i,j});
        end
        if strcmp(modeflag,'classification') && decomposeflag ~= 9
            % First generate label/index vectors for binary classification 
            % for each outer training fold:
            % 1) the training data:
            % label(cv.TrainInd{i,j}) = these are the subject labels of the
            % outer training partion [i,j]
            switch appendfl
                case {0,1,3}
                    
                    if appendfl == 3
                        cvin = cv.cvin{i,j};
                    else
                        cvin = [];
                    end
                    
                    cv.class{i,j} = GenClass(g, ulb, nclass, Label(cv.TrainInd{i,j}), decomposeflag, NaNflag);

                    if ~isempty(Constraint)
                        [cv.class{i,j}, cv.cvin{i,j}, InFold] = makefolds(Label(cv.TrainInd{i,j}), ...
                            cv.class{i,j}, InFold, InPerms, decomposeflag, Constraint(cv.TrainInd{i,j}), ...
                            appendfl , cvin , Eq, [], auto_adjust, cv2frame);
                    else
                        if ~isempty(CV1LCO)
                            fprintf('\nGenerate CV1 leave-group out partitions for outer training partition [%g,%g].',i,j)
                            [cv.class{i,j}, cv.cvin{i,j}, InFold] = makefolds(Label(cv.TrainInd{i,j}), ...
                                cv.class{i,j}, InFold, InPerms, decomposeflag, [], appendfl, cvin, Eq, CV1LCO.ind(cv.TrainInd{i,j}), auto_adjust, cv2frame);
                        else
                            fprintf('\nGenerate CV1 crossvalidation partitions for outer training partition [%g,%g].',i,j)
                            [cv.class{i,j}, cv.cvin{i,j}, InFold] = makefolds(Label(cv.TrainInd{i,j}), ...
                                cv.class{i,j}, InFold, InPerms, decomposeflag, [], appendfl, cvin, Eq, [], auto_adjust, cv2frame);
                        end
                    end
                    % 2) the validation data:
                    % label(cv.TestInd{i,j}) = these are the subject labels of the
                    % outer validation partion [i,j]
                    cv.classnew{i,j} = GenClass(g, ulb, nclass, Label(cv.TestInd{i,j}), decomposeflag, NaNflag);
                    
                case 2
                    
                    if ~isempty(Constraint)
                        [cv.class{i,j}, cv.cvin{i,j}, InFold] = makefolds(Label(cv.TrainInd{i,j}), ...
                            cv.class{i,j}, InFold, InPerms, decomposeflag, Constraint(cv.TrainInd{i,j}), ...
                            appendfl, cv.cvin{i,j}, Eq, [], auto_adjust, cv2frame);
                    else
                        if ~isempty(CV1LCO)
                            fprintf('\nGenerate CV1 leave-group out partitions for outer training partition [%g,%g].',i,j)
                            [cv.class{i,j}, cv.cvin{i,j}, InFold] = makefolds(Label(cv.TrainInd{i,j}), ...
                                cv.class{i,j}, InFold, InPerms, decomposeflag, [], appendfl, cv.cvin{i,j}, Eq, CV1LCO.ind(cv.TrainInd{i,j}), auto_adjust, cv2frame);
                        else
                            fprintf('\nGenerate CV1 crossvalidation partitions for outer training partition [%g,%g].',i,j)
                            [cv.class{i,j}, cv.cvin{i,j}, InFold] = makefolds(Label(cv.TrainInd{i,j}), ...
                                cv.class{i,j}, InFold, InPerms, decomposeflag, [], appendfl, cv.cvin{i,j}, Eq, [], auto_adjust, cv2frame);
                        end
                    end
            end
            
        elseif  strcmp(modeflag,'regression') || (strcmp(modeflag,'classification') && decomposeflag == 9)
            fprintf('\nGenerate CV1 crossvalidation partitions for outer training partition [%g,%g].',i,j)
            
            switch appendfl
               
                case {0,1}
                    if ~isempty(Constraint)
                        cv.cvin{i,j} = nk_CVpartition2(InPerms, InFold, ...
                            Label(cv.TrainInd{i,j}), Constraint(cv.TrainInd{i,j}), Eq, auto_adjust);
                    else
                        if CV1LOO
                            cv.cvin{i,j} = nk_LOOpartition(Label(cv.TrainInd{i,j}));
                        elseif ~isempty(CV1LCO)
                            cv.cvin{i,j} = nk_INDEPpartition(CV1LCO.ind(cv.TrainInd{i,j}), Label(cv.TrainInd{i,j}), cv2frame);
                        else
                            cv.cvin{i,j} = nk_CVpartition2(InPerms, InFold, Label(cv.TrainInd{i,j}), [], Eq, auto_adjust);
                        end
                    end
                case 2
                    if ~isempty(Constraint)
                        cv.cvin{i,j} = nk_CVpartition2(InPerms, InFold, ...
                            Label(cv.TrainInd{i,j}), Constraint(cv.TrainInd{i,j}), Eq, auto_adjust);
                    else
                        if CV1LOO
                            cv.cvin{i,j} = nk_LOOpartition(Label(cv.TrainInd{i,j}));
                        elseif ~isempty(CV1LCO)
                            cv.cvin{i,j} = nk_INDEPpartition(CV1LCO.ind(cv.TrainInd{i,j}), Label(cv.TrainInd{i,j}), cv2frame);
                        else
                            cv.cvin{i,j} = nk_CVpartition2(InPerms, InFold, Label(cv.TrainInd{i,j}), [], Eq, auto_adjust);
                        end
                    end
            end
        end
    end
end

return
% ______________________________
function class = GenClass(g, xlb, nclass, lb, decomposeflag, nanflag)

switch decomposeflag
    
    % Generate binary classification classes (One-Vs-One)
    case 1
        
        cnt=1; class = cell(nclass*(nclass-1)/2,1);
        
        for i=1:nclass-1
            
            ijlb = zeros(1,2);
            ijlb(1) = xlb(i);
            
            for j=i+1:nclass
                
                ijlb(2) = xlb(j);
                class{cnt}.groups = ijlb;
                
                if ~isempty(g{i} ) && ~isempty(g{j})
                    class{cnt}.groupdesc = [g{i} ' vs ' g{j}];
                end
                
                % Positive label
                ind1 = find( lb == ijlb(1) );
                ind2 = find( lb == ijlb(2) );
                label1 = ones(1,numel(ind1))';
                label2 = -1*ones(1,numel(ind2))';
                class{cnt}.ind = [ind1; ind2];
                class{cnt}.label = [ label1; label2 ];
                if nanflag,
                    indnan = find( ~isfinite(lb) );
                    labelnan = nan(numel(indnan),1);
                    class{cnt}.ind = [class{cnt}.ind; indnan];
                    class{cnt}.label = [ class{cnt}.label; labelnan ];
                end
                cnt=cnt+1;
            end
        end
    
    % Generate binary classification classes (One-Vs-All)
    case 2
        class = cell(nclass,1);
        for i=1:nclass
            if ~isempty(g{i}), class{i}.groupdesc = [g{i} ' vs ALL']; end
            class{i}.groups(1) = xlb(i);
            indpos = lb==xlb(i); indneg = lb~=xlb(i);
            label = zeros(size(lb)); 
            if any(indpos),label(indpos) = 1; end
            if any(indneg),label(indneg) = -1; end
            class{i}.label = label;
            class{i}.ind = (1:size(lb,1))';
            if nanflag,
                indnan = find( ~isfinite(lb) );
                labelnan = nan(numel(indnan),1);
                class{i}.ind = [class{i}.ind; indnan];
                class{i}.label = [ class{i}.label; labelnan ];
            end
        end
    
    % Generate multi-group classification setup
    case 9
        class.groups =  1 : nclass ;
        class.groupdesc = 'Multi-group classification';
        class.label = lb;
        class.ind = (1:size(lb,1))';
        if nanflag,
            indnan = find( ~isfinite(lb) );
            labelnan = nan(numel(indnan),1);
            class.ind = [class.ind; indnan];
            class.label = [ class.label; labelnan ];
        end
end
	
% ___________________________________________
function [class, cv, xfolds] = makefolds(label, class, nfolds, nperms, decomposeflag, constraint,  appendfl, oldcv, eq, lco, auto_adjust, cv2frame)

groupind = [1,-1]; if ~exist('auto_adjust', 'var'), auto_adjust = false; end
xfolds = nfolds;
% Generate CV1 partitions
if appendfl ~=3
    if exist('lco','var') && ~isempty(lco)
        [cv, nfolds, nperms] = nk_INDEPpartition(lco, label, cv2frame);
    else
        LOOflag = false; if nfolds == -1; LOOflag = true; end
        if ~LOOflag
            cv = nk_CVpartition2(nperms, nfolds, label, constraint, eq, auto_adjust);
            if ~isstruct(cv);
                nfolds = cv; cv = nk_CVpartition2(nperms, nfolds, label, constraint, eq, auto_adjust);
            end
        else
            nfolds = numel(label); nperms = 1;
            cv = nk_LOOpartition(label);
        end
    end
else
    cv = oldcv; clear oldcv;
end
tclass = class;

for i=1:length(class)
    
    if isfield(tclass,'TrainInd')
        tclass{i}.TrainInd{nperms,nfolds}   = [];
        tclass{i}.TestInd{nperms,nfolds}    = [];
        tclass{i}.TrainLabel{nperms,nfolds} = [];
        tclass{i}.TestLabel{nperms,nfolds}  = [];
        %tclass{i}.TestFoldInd{nperms}=[];
    else
        if iscell(tclass)
            tclass{i}.TrainInd                  = cell(nperms,nfolds);
            tclass{i}.TestInd                   = cell(nperms,nfolds);
            tclass{i}.TrainLabel                = cell(nperms,nfolds);
            tclass{i}.TestLabel                 = cell(nperms,nfolds);
        else
            tclass.TrainInd                  = cell(nperms,nfolds);
            tclass.TestInd                   = cell(nperms,nfolds);
            tclass.TrainLabel                = cell(nperms,nfolds);
            tclass.TestLabel                 = cell(nperms,nfolds);
        end
    end    
    
    for j=1:nperms

        %tclass{i}.TestFoldInd{j}=zeros(size(class{i}.ind));

        for k=1:nfolds
 
                switch decomposeflag
                    
                    case 1 % One-vs-One
                        
                        for l=1:length(groupind)
                        
                            %tclass{i}.TestFoldInd{j}(cv.TestInd{j,k}(label(cv.TestInd{j,k})
                            %== class{i}.groups(l)))=k;
                            belongstrain = find(label(cv.TrainInd{j,k}) == class{i}.groups(l));
                            belongstest = find(label(cv.TestInd{j,k}) == class{i}.groups(l));

                            if ~isempty(belongstrain)
                                tclass{i}.TrainInd{j,k} = [tclass{i}.TrainInd{j,k}; ...
                                    cv.TrainInd{j,k}(belongstrain)];
                                tclass{i}.TrainLabel{j,k} = [tclass{i}.TrainLabel{j,k}; ...
                                    groupind(l)*ones(size(belongstrain,1),1)];
                            end

                            if ~isempty(belongstest)
                                tclass{i}.TestInd{j,k} = [tclass{i}.TestInd{j,k}; ...
                                     cv.TestInd{j,k}(belongstest)];
                                tclass{i}.TestLabel{j,k} = [tclass{i}.TestLabel{j,k}; ...
                                    groupind(l)*ones(size(belongstest,1),1)];
                            end
                            
                        end
                        belongstrain = find(~isfinite(label(cv.TrainInd{j,k})));
                        if ~isempty(belongstrain)
                            tclass{i}.TrainInd{j,k} = [tclass{i}.TrainInd{j,k}; ...
                                cv.TrainInd{j,k}(belongstrain)];
                            tclass{i}.TrainLabel{j,k} = [tclass{i}.TrainLabel{j,k}; ...
                                nan(size(belongstrain,1),1)];
                        end
                        
                    case 2 % One-Vs-All
                        
                        %tclass{i}.TestFoldInd{j}(cv.TestInd{j,k}(label(cv.TestInd{j,k}) == class{i}.groups))=k;
                        tclass{i}.TrainInd{j,k}  = cv.TrainInd{j,k};
                        lb = zeros(size(cv.TrainInd{j,k}));
                        indpos = label(cv.TrainInd{j,k}) == class{i}.groups; 
                        indneg = label(cv.TrainInd{j,k}) ~= class{i}.groups;
                        indnan = ~isfinite(label(cv.TrainInd{j,k}));
                        lb(indpos) = 1; lb(indneg) = -1;
                        if ~isempty(indnan), lb(indnan) = NaN; end
                        tclass{i}.TrainLabel{j,k} = lb;
                        
                        tclass{i}.TestInd{j,k}   = cv.TestInd{j,k};
                        lb = zeros(size(cv.TestInd{j,k}));
                        indpos = label(cv.TestInd{j,k}) == class{i}.groups; 
                        indneg = label(cv.TestInd{j,k}) ~= class{i}.groups;
                        lb(indpos) = 1; lb(indneg) = -1;
                        tclass{i}.TestLabel{j,k} = lb;
                        
                    case 9 % Multi-group

                        tclass.TrainInd{j,k}  = cv.TrainInd{j,k};
                        tclass.TestInd{j,k}   = cv.TestInd{j,k};
                        tclass.TrainLabel{j,k} = label(cv.TrainInd{j,k});
                        tclass.TestLabel{j,k} = label(cv.TestInd{j,k});
                end
        end
        
        %tclass{i}.TestFoldInd{j} =  uint8(tclass{i}.TestFoldInd{j}(tclass{i}.TestFoldInd{j}~=0));
    end
    if iscell(tclass)
        switch appendfl
            case {0,1,3}
                class{i} = tclass{i};
            case 2
                class{i}.TrainInd = [ class{i}.TrainInd; tclass{i}.TrainInd ];
                class{i}.TestInd = [ class{i}.TestInd; tclass{i}.TestInd ];
                class{i}.TrainLabel = [ class{i}.TrainLabel; tclass{i}.TrainLabel ];
                class{i}.TestLabel = [ class{i}.TestLabel; tclass{i}.TestLabel ];
                %class{i}.TestFoldInd = [ class{i}.TestFoldInd tclass{i}.TestFoldInd ];
        end
    else
        switch appendfl
            case {0,1,3}
                class = tclass;
            case 2
                class.TrainInd = [ class.TrainInd; tclass.TrainInd ];
                class.TestInd = [ class.TestInd; tclass.TestInd ];
                class.TrainLabel = [ class.TrainLabel; tclass.TrainLabel ];
                class.TestLabel = [ class.TestLabel; tclass.TestLabel ];
                %class.TestFoldInd = [ class.TestFoldInd tclass.TestFoldInd ];
        end
    end
end

if exist('oldcv','var') && ~isempty(oldcv) && appendfl == 2
     cv.TrainInd = [ oldcv.TrainInd; cv.TrainInd];
     cv.TestInd = [ oldcv.TestInd; cv.TestInd];
end

return
