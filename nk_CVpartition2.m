% =========================================================================
% FORMAT cv = nk_CVpartition2(nperms, K, Labels, Constraint)
% =========================================================================
% This function generates nperms*K partitions of data
% using cross-validation resampling
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Nikolaos Koutsouleris 01/2012
function cv = nk_CVpartition2(nperms, K, Labels, Constraint, Eq, AutoAdjust)
global NM

trainidxs = cell(nperms,K); testidxs = cell(nperms,K); 

uLabels = unique(Labels);
% Check whether there are non-finite values in the label vector
if any(~isfinite(uLabels)), 
    NaNflag = true; 
    uLabels(~isfinite(uLabels))=[];
else
    NaNflag = false;
end

mL = numel(uLabels);
if ~exist('Eq','var'), Eq = []; end
if ~exist('Constraint','var') || isempty(Constraint)
    mC = 1; uConstraint = []; Constraint = [];
else
    uConstraint = unique(Constraint);
    mC = numel(uConstraint); 
end

if ~isempty(uConstraint)
    C = zeros(mL,mC);
    for j = 1:mL
       for hu = 1:mC
           C(j,hu) = sum( Labels == uLabels(j) & Constraint == uConstraint(hu) );
       end
    end
    minC = min(C,[],2);
end  

% Generate Permutation indices
permmat = nk_PermInd2(nperms, Labels, Constraint);

for h=1:nperms % Loop through perms
    
    rInd        = permmat(h,:)';
    trainidx    = cell(1,K);
    testidx     = cell(1,K);
    
    for j=1:mL,  % Loop through classes
        
        % row vector, indeces of members of the current (jth) class
        indClassCX = []; 
        indLabels = find(Labels == uLabels(j));
        
        if ~isempty(uConstraint)
            indClass = [];
            for hu=1:mC
                indClassX = find(Labels == uLabels(j) & Constraint == uConstraint(hu));
                indClass =  [indClass; indClassX( 1 : minC(j)) ];
            end
            indRem = setdiff(indLabels,indClass);
            ConstrXClass = Constraint(indClass);
            if (numel(ConstrXClass) / K) < numel(uConstraint)
                % There are less subjects in current constraint x class group then K folds
                CXfl = true;
                fprintf('\tNot enough observations in current constraint x group class')
            else
                CXfl = false;
            end
        else
            indClass = indLabels;
            indRem = [];
        end
        
        nClassMem = length(indClass);
        testsize = floor(nClassMem/K);

        for i=1:K % Loop through folds
             
            % number of members (examples) in the class    
            % ie: nmem >= K
            if testsize > 0,
                if ~isempty(uConstraint)
                   
                    % Loop through constraint classes
                    endpos = 1;
                    for hu = 1:mC
                        indC = find(ConstrXClass == uConstraint(hu));
                        if ~CXfl, 
                            testsizeC = floor(numel(indC)/K);
                        else
                            testsizeC = 1;
                        end
                        startpos = (i-1)*testsizeC + 1;
                        endpos = i*testsizeC;
                        if CXfl
                            fprintf('\tRepeating test subjects at fold %g', i)
                            if endpos > numel(indC); endpos = numel(indC); end
                            if startpos > numel(indC); startpos = numel(indC); end
                        end
                        indCx = indC(startpos:endpos);
                        testidx{i} = [testidx{i}; indClass(indCx)];
                        indClassCX = [indClassCX; indClass(indCx)];
                        
                    end
                else
                    startpos = (i-1)*testsize +1; 
                    endpos = i*testsize;
                    testidx{i} = [testidx{i}; indClass(startpos:endpos)];
                    indClassCX = [indClassCX; indClass(startpos:endpos)];
                    % In case of regression equalize target label histogram
                    % by removing subjects from overrepresented label bins
                end
                
            else
                if exist('AutoAdjust','var') && AutoAdjust 
                    cv = nClassMem; return;
                else
                    AdjStr = ['Adjust to N<=' num2str(nClassMem)];
                    choice = questdlg(sprintf('Not enough members to send to test data.\nNumber of CV2 folds has to be <=%g',nClassMem), ...
                        'Cross-validation structure generator','Abort', AdjStr,2);
                    switch choice
                        case 'Abort'
                            error('Cross-validation generation aborted. Adjust your settings!')
                        case AdjStr
                            cv = nClassMem; return;
                    end
                end
                % We let all points of this class participate in training but no point in test:
                % We could have let only the ith member of this class (if any) be for test and the rest for train.
            end
            end
        
        %% Create stratified cross-validation scheme
        if ~isempty(indClassCX) 
            
            indRem = [indRem; setdiff(indClass,indClassCX)];
            cnt = numel(indRem); pInd = []; 
            if cnt > K
                cntK = K;
            else
                cntK = cnt;
            end
            while cnt > 0
                indX = randperm( K );
                pInd = [pInd indX(1:cntK)];
                cnt = cnt - cntK;
                if cnt <= K
                    cntK = cnt;
                end
            end
            for i = 1:numel(indRem)
                testidx{pInd(i)} = [testidx{pInd(i)}; indRem(i)];
            end
        end
        for i = 1:K
            trainidx{i} = [trainidx{i}; setdiff(indLabels, testidx{i})];
        end        
    end
    % Add NaN labels to each training population
    if NaNflag
        for i = 1:K
            trainidx{i} = [trainidx{i}; find(~isfinite(Labels))];
        end        
    end
    %% Do balancing / histogram equalization in each fold
    if ~isempty(Eq)
        for i=1:K
            [ removed, retained ] = nk_EqualizeHisto(Eq, Eq.Covar(trainidx{i}), trainidx{i}, NM.modeflag);
            trainidx{i} = retained;
            if Eq.AddRemoved2Test, testidx{i} = [testidx{i}; removed]; end
        end
    end           
    % Convert to integer to save space
    for i=1:K
        testidxs{h,i} = uint16(rInd(testidx{i}));
        trainidxs{h,i} = uint16(rInd(trainidx{i}));
    end
end

cv.TrainInd = trainidxs;
cv.TestInd = testidxs;

return
