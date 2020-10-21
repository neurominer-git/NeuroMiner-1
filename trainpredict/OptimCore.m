% =========================================================================
% FORMAT Res = OptimCore(IN, OUT, strout, Param, OptMode, Weighting)
% =========================================================================
% The OPTIMCORE function is the central pacemaker of the training module of 
% NM. It allows to invoke different optimization strategies such as filters 
% or wrappers (via FoldPerm). This includes ensemble-based optimization 
% strategies tailored to binary/regression and multi-class settings. 
% Furthermore, it collects the CV1 training and CV1 test results 
% across the CV1 partitions and organizes them into a structured output.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%(c) Nikolaos Koutsouleris, 01/2019

function Res = OptimCore(IN, OUT, strout, Param, OptMode, Weighting)

global SVM MULTI MODEFL RFE VERBOSE
                        
% STEP 1a)
% Train prediction model on pre-filtered features subspaces with an
% optional wrapper-type optimization (OptMode = 1) of these feature 
% subspaces:
% CV1-Training, CV1-Test data performance (Tr, Ts) are returned as well
% as the predicted target values / labels for the training and test
% data. In linear SVMs also margin params are returned.

rtrfl = 0;
if isfield(RFE,'RetrainImmediate') && ...
        RFE.RetrainImmediate && ...
        RFE.ClassRetrain && ...
        (~RFE.Wrapper.flag || OptMode == 1)
    rtrfl = 1;
end

if rtrfl
    % Retrain without optimization
    [ IN, OUT ] = FoldPerm(IN, OUT, 'Immediate Retrain (CV1-Tr + CV1-Ts)', 0, 1, 1, 0);
else
    % Train models with optimization
    [ IN, OUT ] = FoldPerm(IN, OUT, strout, OptMode, 0, 0, Param.SubSpaceStepping);
end

% STEP 1b)
% If Multi-group optimization has to be carried out, perform multi-group
% feature subset selection using obtained decision values /
% probabilities
if MULTI.flag
    strprf = 'Multi-Group classifier'; OUT = MultiDichoFoldPerm(IN, OUT);
else
    switch MODEFL
        case 'classification'
            strprf = 'Binary classifier';
        case 'regression'
            strprf = 'Regressor';
    end
end

% STEP 2)
% Select subspaces according to subspace performance 
% --------------------------------------------------
% in each data partition following to different strategies
% 1) "winner takes it all"
% 2) "within range from max"
% 3) "above percentile"

if VERBOSE
    fprintf('\n\nEvaluate CV performance and select subsets')
    fprintf('\n------------------------------------------')
end
strout = [strprf ' selection'];
minmaxfl = 1; 
RetrainFlag = RFE.ClassRetrain && (~RFE.Wrapper.flag || OptMode == 1) && ~rtrfl;
switch Param.CostFun
    case 0
        tRankCrit = OUT.tr; 
        fprintf('\nNo subspace evaluation required.')

    case {1,2,3} % Maximize the performance of SVM on the CV1 test data
        switch Param.CostFun
            case 1
                 tRankCrit = OUT.tr;
                 if MULTI.flag && MULTI.train, tRankCrit = OUT.mtr; end
            case 2
                 tRankCrit = OUT.ts;
                 if MULTI.flag && MULTI.train, tRankCrit = OUT.mts; end
            case 3
                 if MULTI.flag && MULTI.train, 
                     X = cellfun(@plus,OUT.mts,OUT.mtr,'UniformOutput',false);
                     tRankCrit = cellfun(@rdivide,X, repmat({2},size(X)),'UniformOutput',false);
                 else
                     X = cellfun(@plus,OUT.ts,OUT.tr,'UniformOutput',false);
                     tRankCrit = cellfun(@rdivide,X, repmat({2},size(X)),'UniformOutput',false);
                 end
        end
        if ~strcmp(MODEFL,'classification') && any(SVM.GridParam == [ 9 11 12 18 ]), minmaxfl = 2; end
end

% Evaluate margin based criterion / performance in all feature subspaces
% and extract one (or more) subspaces that fulfill the criterion defined
% in SubSpaceStrategy (see there)

% Optimize using binary performance
[OUT.F, OUT.Weights] = EvalSubSpaces(tRankCrit, strout, ...
          Param.SubSpaceStrategy, Param.SubSpaceCrit, minmaxfl, [], Weighting);

% Optionally, detrend predictions errors in case of regression
OUT = DetrendFoldPerm(IN, OUT);

% if some sort of ensemble predictor should be constructed, do that
% at this point using definitions in RFE.Filter.EnsembleStrategy
strout='Construct';
[ IN , OUT ] = SubSpaceStrat(IN, OUT, Param, OptMode, strout);

% Cross-CV1 probabilistic feature selection
if isfield(Param,'PFE') && Param.PFE.flag
    [ IN , OUT ] = SelectFeaturesAcrossCV1(IN, OUT, Param, minmaxfl, 0);
end

% Optionally retrain the classifier by conactenating CV1 training and test data
S=zeros(IN.nclass,1);

if RetrainFlag
     [ IN, OUT ] = FoldPerm(IN, OUT, 'Retrain (CV1-Tr + CV1-Ts)', 0, 1, 0, 0);
     for curclass=1:IN.nclass
        S(curclass) = size(IN.Y.TrL{1,1}{curclass},1) + size(IN.Y.CVL{1,1}{curclass},1);
     end
elseif rtrfl
    for curclass=1:IN.nclass
        S(curclass) = size(IN.Y.TrL{1,1}{curclass},1) + size(IN.Y.CVL{1,1}{curclass},1);
    end
else
    for curclass=1:IN.nclass
        S(curclass) = size(IN.Y.TrL{1,1}{curclass},1);
    end
end

% Evaluate Sequence Optimizer
switch SVM.prog
    case 'SEQOPT'
        [OUT.critgain, OUT.examfreq, OUT.percthreshU, OUT.percthreshL, OUT.absthreshU, OUT.absthreshL] = EvalSeqOpt(OUT.mdl);    
    case 'WBLCOX'
        [OUT.threshperc, OUT.threshprob, OUT.times] = EvalCoxPH(OUT.mdl);
end

% Copy results to output structure
Res = TransferResults(IN, OUT, Param);

% Compute model complexity parameters
[Res.ModelComplexity.sumNSV, ...
Res.ModelComplexity.meanNSV, ...
Res.ModelComplexity.stdNSV, ...
Res.ModelComplexity.Complex] = nk_ModelComplexityParams(Res.Models, S);

end
