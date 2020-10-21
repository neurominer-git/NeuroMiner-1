function r = rfe_algo_settings_multi(Y, mY, label, labelB, labelM, mYnew, labelnew, labelnewM, Ps, FullFeat, FullParam, ngroups, ActStr)

global RFE VERBOSE SVM TRAINFUNC

r.nclass = numel(Y);
r.ngroups = ngroups;
r.FeatRandPerc = 0;
r.FeatStepPerc = true;

if isfield(RFE.Wrapper.GreedySearch,'PreSort') && RFE.Wrapper.GreedySearch.PreSort, r.PreSort= RFE.Wrapper.GreedySearch.PreSort; else r.PreSort = false; end

for curclass = 1:r.nclass
    r.FullInd{curclass} = ( find(any(Y{curclass}) & std(Y{curclass})~=0 & sum(isfinite(Y{curclass}))~=0 ))' ;
    r.Y{curclass} = Y{curclass}(:,r.FullInd{curclass});
    r.mY{curclass} = mY{curclass}(:,r.FullInd{curclass});
    r.mYnew{curclass} = mYnew{curclass}(:,r.FullInd{curclass}); 
    r.kFea(curclass) = numel(r.FullInd{curclass});
    r.BinaryLabel{curclass} = labelB{curclass}(labelB{curclass}~=0);
    [~, r.FullModel{curclass}] = feval(TRAINFUNC, r.Y{curclass}, r.BinaryLabel{curclass}, 1, Ps{curclass});    
    if r.PreSort
        %r.BinaryDummy{curclass} = single(nk_MakeDummyVariables(r.BinaryLabel{curclass}));
        [~, r.FullModel{curclass}] = feval(TRAINFUNC, r.Y{curclass}, r.BinaryLabel{curclass}, 1, Ps{curclass});    
        r.W{curclass} = nk_GetAlgoWeightVec(SVM, r.Y{curclass}, r.BinaryLabel{curclass}, r.FullModel{curclass}, true, 0);
        if isempty(r.W{curclass}), r.W{curclass} = nk_CorrMat(r.Y{curclass},r.BinaryLabel{curclass}); end
        [~,r.PreOrder{curclass}] = sort(abs(r.W{curclass}),'ascend');
    end   
end

switch RFE.Wrapper.type
    case 1
        %% Feature Sorting Mode
        % if WeightSort==2 then the feature weight vector will be used for sorting
        % the features. This is only valid for linear models

        %% Feature block size settings
        for curclass=1:numel(Y)
            kFea = length(r.FullInd{curclass}); 
            if ~isfield(RFE.Wrapper.GreedySearch,'FeatStepPerc')
                r.lperc(curclass) = 1/(kFea/100); % lperc should multiply to a feature stepping of 1 
                r.FeatStepPerc(curclass) = false;
            elseif ~RFE.Wrapper.GreedySearch.FeatStepPerc
                r.lperc(curclass) = 1;
            else
                r.lperc(curclass) = RFE.Wrapper.GreedySearch.FeatStepPerc;
            end

            %% Early stopping setting
            r.MinNum(curclass) = 1;
            if RFE.Wrapper.GreedySearch.EarlyStop.Thresh
               switch RFE.Wrapper.GreedySearch.EarlyStop.Perc
                   case 1 % Percentage Mode
                       r.MinNum(curclass) = size(r.Y{curclass},2) / 100 * RFE.Wrapper.GreedySearch.EarlyStop.Thresh; 
                   case 2 % Absolute mMde
                       r.MinNum(curclass) = RFE.Wrapper.GreedySearch.EarlyStop.Thresh;
               end
            end
        end
        if isfield(RFE.Wrapper.GreedySearch,'FeatRandPerc')
            r.FeatRandPerc = RFE.Wrapper.GreedySearch.FeatRandPerc;
        end
        r.perm = false;
        if isfield(RFE.Wrapper.GreedySearch,'PERM') && RFE.Wrapper.GreedySearch.PERM.flag
            r.perm = true;
            r.nperms = 100;
        end
end

% Optimization criterion setup
[r.evaldir, ~, r.optfunc, r.optparam] = nk_ReturnEvalOperator(SVM.GridParam);

%% Obtain full feature space performance, if not already done so
if ~exist('FullParam','var') || isempty(FullParam)
    
    switch ActStr
        case 'Tr'
            r.T = r.Y;
            %r.mT = r.mY;
            r.L = label;
            r.Lm = labelM;
        case 'CV'
            r.T = r.mYnew;
            %r.mT = r.mYnew;
            r.L = labelnew;
            r.Lm = labelnewM;
        case 'TrCV'
            for curclass=1:numel(Y)
                r.T{curclass} = [r.mY{curclass}; r.mYnew{curclass}];
                r.L{curclass} = [labelB{curclass}; labelnew{curclass}];
            end
            r.Lm = [labelM; labelnewM];

    end
    ds = zeros(size(r.T{1},1),1);
    ts = zeros(size(r.T{1},1),1);
     for curclass=1:numel(Y)
         if isempty(r.FullModel{curclass})
            [~, r.FullModel{curclass}] = feval(TRAINFUNC, r.Y{curclass}, label{curclass}, 1, Ps{curclass});    
         end
         [r.FullParam{curclass}, ds(:,curclass), ts(:,curclass) ]= nk_GetTestPerf(r.T{curclass}, r.L{curclass}, [], r.FullModel{curclass}, r.Y{curclass});
         if VERBOSE, fprintf('\nFull model #%g:\t# Features: %g, %s = %g', ...
                            curclass, numel(r.FullInd{curclass}), ActStr, r.FullParam{curclass}), end
     end
     r.FullParamMulti = nk_MultiEnsPerf(ds, ts, r.Lm, 1:numel(Y), r.ngroups);
     if VERBOSE, fprintf('\nFull multi-class model: %s = %g', ActStr, r.FullParamMulti); end
     
else
    r.FullParam = FullParam;
    r.FullModel = [];
end

if isfield(RFE.Wrapper.GreedySearch,'KneePointDetection') && RFE.Wrapper.GreedySearch.KneePointDetection == 1
    r.KneePoint = true;
else
    r.KneePoint = false;
end

if isfield(RFE.Wrapper.GreedySearch,'MultiClassOptimization') && RFE.Wrapper.GreedySearch.MultiClassOptimization == 1
    if RFE.Wrapper.GreedySearch.PERM.flag
        r.perm = true;
        r.nperms = RFE.Wrapper.GreedySearch.PERM.nperms;
    end
end

if isfield(RFE.Wrapper.GreedySearch,'CaseFrwd'),
    r.CaseU = RFE.Wrapper.GreedySearch.CaseFrwd.Upper;
    r.CaseL = RFE.Wrapper.GreedySearch.CaseFrwd.Lower;
    r.CaseStep = RFE.Wrapper.GreedySearch.CaseFrwd.Step;
    r.CaseSpU = r.CaseU(2):-1*round((r.CaseU(2)-r.CaseU(1))/r.CaseStep):r.CaseU(1);
    r.CaseSpL = r.CaseL(1):round((r.CaseL(2)-r.CaseL(1))/r.CaseStep):r.CaseL(2);
end
