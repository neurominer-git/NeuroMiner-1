%==========================================================================
%FORMAT [ IN , OUT ] = SubSpaceStrat(IN, OUT, Param.EnsembleStrategy, strout)
%==========================================================================
%Performs subspace optimization using either probabilistic feature      
%extraction or ensemble-based optimization                             
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%(c) Nikolaos Koutsouleris, 04/2015

function [ IN , OUT ] = SubSpaceStrat(IN, OUT, Param, OptMode, strout)
                  
global CV MULTI MODEFL W2AVAIL MULTILABEL VERBOSE

%%%% INITIALIZE %%%%

% Prepare cell containers for ensemble learning
OUT.TrHDperf    = zeros( IN.nperms, IN.nfolds, IN.nclass );
OUT.TrHTperf    = zeros( IN.nperms, IN.nfolds, IN.nclass );
OUT.CVHDperf    = zeros( IN.nperms, IN.nfolds, IN.nclass );
OUT.CVHTperf    = zeros( IN.nperms, IN.nfolds, IN.nclass );
OUT.TrDiv       = zeros( IN.nperms, IN.nfolds, IN.nclass );
OUT.CVDiv       = zeros( IN.nperms, IN.nfolds, IN.nclass );
OUT.TrHD        = cell( IN.nperms, IN.nfolds, IN.nclass );
OUT.TrHT        = cell( IN.nperms, IN.nfolds, IN.nclass );
OUT.CVHD        = cell( IN.nperms, IN.nfolds, IN.nclass );
OUT.CVHT        = cell( IN.nperms, IN.nfolds, IN.nclass );

nc = size(OUT.F,3);
ns = size(IN.F,3);
EnsType = 0;
Metric = 2;
DivFunc = 'nk_Entropy';
ngroups = IN.ngroups;
if isfield(Param,'EnsembleStrategy') && Param.SubSpaceFlag 
    if ~isempty(Param.EnsembleStrategy) && isfield(Param.EnsembleStrategy,'type')
        EnsType = Param.EnsembleStrategy.type;
        Metric = Param.EnsembleStrategy.Metric;
        DivFunc = Param.EnsembleStrategy.DivFunc;
    end
end

if MULTI.flag
    % If multi-group training option has been chosen prepare multi-group
    % prediction and peformance copntainers
    OUT.mTrPred     = cell( IN.nperms, IN.nfolds );
    OUT.mCVPred     = cell( IN.nperms, IN.nfolds );
    OUT.mTrPerf     = zeros( IN.nperms, IN.nfolds );
    OUT.mCVPerf     = zeros( IN.nperms, IN.nfolds );
    OUT.mTrDiv      = zeros( IN.nperms, IN.nfolds );
    OUT.mCVDiv      = zeros( IN.nperms, IN.nfolds );
    tF        = cell( IN.nperms, IN.nfolds, IN.nclass );
    tWeights  = cell( IN.nperms, IN.nfolds, IN.nclass );
    tSubSpaces  = cell( IN.nperms, IN.nfolds, ns );
   
else
    tF          = cell( IN.nperms, IN.nfolds, IN.nclass );
    tWeights    = cell( IN.nperms, IN.nfolds, IN.nclass );
    tSubSpaces  = cell( IN.nperms, IN.nfolds, IN.nclass );
end

tModels         = cell( IN.nperms, IN.nfolds, IN.nclass );

if W2AVAIL
    tW2         = cell( IN.nperms, IN.nfolds, IN.nclass );
    tMd         = cell( IN.nperms, IN.nfolds, IN.nclass );
    tMm         = cell( IN.nperms, IN.nfolds, IN.nclass );
end

%%%% APPLY SUBSPACE STRATEGY %%%%

for i=1:IN.nperms % loop through CV1 permutations
    
    for j=1:IN.nfolds  % loop through CV1 folds
              
        mTr = []; mCV = []; mC = []; tkIndCat = [];
        kIndN = cell(IN.nclass,1); TrL = kIndN; CVL = kIndN;

        switch EnsType

            case 9 

                %% PROBABILISTIC FEATURE SPACE CONSTRUCTION

                for curclass=1:IN.nclass

                    %% Get feature subspace mask and feature weighting:
                    if nc > 1 
                        % Every dichotomizer has its own feature subspace mask
                        kInd = OUT.F{i,j,curclass}; W = OUT.Weights{i,j,curclass};
                    else
                        % All dichotomizers share one feature subspace mask
                        % This mode is used for multi-group optimization
                        kInd = OUT.F{i,j}; W = OUT.Weights{i,j};
                        ndim = size(OUT.Trdecs{i,j,curclass},2);
                         % number of base learners is not equal to the number of features in feature mask
                        if numel(kInd) > ndim, kInd = true(ndim,1); end 
                    end

                    % => only one classifier/predictor will be trained
                    if nc == 1 
                        OUT.F{i,j}       = 1;
                        OUT.Weights{i,j} = 1;
                        if ns > 1
                            IN.F{i,j,curclass} = ProbabilisticFea( IN.F{i,j,curclass}(:,OUT.kxVec{i,j,curclass}(kInd)), Param.EnsembleStrategy );
                        else
                            IN.F{i,j} = ProbabilisticFea( IN.F{i,j}(:,OUT.kxVec{i,j,curclass}(kInd)), Param.EnsembleStrategy );
                        end
                    else
                        OUT.F{i,j,curclass} = 1;
                        OUT.Weights{i,j,curclass} = 1;
                        IN.F{i,j,curclass} = ProbabilisticFea( IN.F{i,j,curclass}(:,OUT.kxVec{i,j,curclass}(kInd)), Param.EnsembleStrategy );
                    end    
                    % Retrain using new feature subspace mask
                    LoopParam.PermVec = i; LoopParam.FoldVec = j; LoopParam.ClassVec = curclass;
                    [IN, OUT] = FoldPerm(IN, OUT, 'Retrain for PFC', OptMode, 0, 0, Param.SubSpaceStepping, LoopParam); 

                    OUT.TrHDperf(i,j,curclass) = mean(OUT.tr{i,j,curclass});
                    OUT.TrHTperf(i,j,curclass) = mean(OUT.tr{i,j,curclass});       
                    OUT.CVHDperf(i,j,curclass) = mean(OUT.ts{i,j,curclass});
                    OUT.CVHTperf(i,j,curclass) = mean(OUT.ts{i,j,curclass}); 

                end

                %% Extract decision values / probabilities using feature subspace mask
                OUT.TrHD(i,j,:) = OUT.Trdecs(i,j,:); 
                OUT.TrHT(i,j,:) = OUT.Trtargs(i,j,:);
                OUT.CVHD(i,j,:) = OUT.CVdecs(i,j,:); 
                OUT.CVHT(i,j,:) = OUT.CVtargs(i,j,:);

                if MULTI.flag,
                    OUT = MultiDichoFoldPerm(IN, OUT, [], LoopParam);
                    OUT.mTrPerf(i,j) = OUT.mtr{i,j};
                    OUT.mCVPerf(i,j) = OUT.mts{i,j};
                end
                
            otherwise

                for curclass = 1: IN.nclass

                    %% Get pointers and labels
                    TrInd                   = IN.Y.TrInd{i,j}{curclass};
                    CVInd                   = IN.Y.CVInd{i,j}{curclass};

                    %% Construct dichotomization labels for CV1 training and test data:
                    CVL{curclass}           = zeros(size(OUT.CVtargs{i,j,curclass},1),1); 
                    CVL{curclass}(CVInd)    = IN.Y.CVL{i,j}{curclass}(:,MULTILABEL.curdim);
                    TrL{curclass}           = zeros(size(OUT.Trtargs{i,j,curclass},1),1); 
                    TrL{curclass}(TrInd)    = IN.Y.TrL{i,j}{curclass}(:,MULTILABEL.curdim);

                    %% Get feature subspace mask and feature weighting:
                    if nc > 1 
                        % Every dichotomizer has its own feature subspace mask
                        kInd = OUT.F{i,j,curclass}; W{curclass} = OUT.Weights{i,j,curclass};
                    else
                        % All dichotomizers share one feature subspace mask
                        % This mode is used for multi-group optimization
                        kInd = OUT.F{i,j}; W{curclass} = OUT.Weights{i,j};
                        ndim = size(OUT.Trdecs{i,j,curclass},2);
                         % number of base learners is not equal to the number of features in feature mask
                        if numel(kInd) > ndim, kInd = OUT.F{i,j}(1:ndim); end 
                    end

                    % Indices to preselected feature subspaces for each
                    % dichotomizer are retrieved using find function
                    kIndN{curclass} = find(kInd); 

                    %% Extract decision values / probabilities using feature subspace mask
                    OUT.TrHD{i,j,curclass} = OUT.Trdecs{i,j,curclass}(:,kIndN{curclass}); 
                    OUT.TrHT{i,j,curclass} = OUT.Trtargs{i,j,curclass}(:,kIndN{curclass});
                    OUT.CVHD{i,j,curclass} = OUT.CVdecs{i,j,curclass}(:,kIndN{curclass}); 
                    OUT.CVHT{i,j,curclass} = OUT.CVtargs{i,j,curclass}(:,kIndN{curclass});

                    %% Weighting: if no weighting is used targets and decision values are multiplied with ones
                    %kW{curclass} = W{curclass}(kIndN{curclass});
                    kW{curclass} = W{curclass}(kIndN{curclass});
                    m1 = size(OUT.TrHD{i,j,curclass},1); 
                    w1 = repmat(kW{curclass},1,m1)';
                    OUT.TrHD{i,j,curclass} = OUT.TrHD{i,j,curclass}.*w1; 
                    OUT.TrHT{i,j,curclass} = OUT.TrHT{i,j,curclass}.*w1;
                    m2 = size(OUT.CVHD{i,j,curclass},1);
                    w2 = repmat(kW{curclass},1,m2)';
                    OUT.CVHD{i,j,curclass} = OUT.CVHD{i,j,curclass}.*w2; 
                    OUT.CVHT{i,j,curclass} = OUT.CVHT{i,j,curclass}.*w2;

                    %% Multi-group preparation
                    % If multi-group classification is performed, perpare
                    % multi-group decision matrix by concatenating decision
                    % values / probabilities across dichotomizers:
                    % Memory allocation is preferable but not implemented
                    % yet.
                    if MULTI.flag
                        % Class array for dichotomization
                        mC = [mC ones(1,size(OUT.CVHT{i,j,curclass},2))*curclass];

                        % Choose either to use decision values / 
                        % probabilities (metric=2) or predicted class memberships 
                        % (metric=1)
                        switch Metric
                            case 1
                                % Training Data
                                mTr = [mTr OUT.TrHT{i,j,curclass}];
                                % Test Data
                                mCV = [mCV OUT.CVHT{i,j,curclass}];
                            case 2
                                mTr = [mTr OUT.TrHD{i,j,curclass}];
                                mCV = [mCV OUT.CVHD{i,j,curclass}];
                        end
                    end

                    switch EnsType

                        case 0 
                            %% AGGREGATED ENSEMBLE (no classifier selection)
                                if size(OUT.TrHD{i,j,curclass},2) > 1 
                                    if VERBOSE, 
                                        switch MODEFL
                                            case 'classification'

                                                fprintf('\n%s (Aggregated Ensemble) => CV1 [%g, %g, %s]:\t', ...
                                                    strout,i,j, CV.class{1,1}{curclass}.groupdesc); 
                                            case 'regression'
                                                fprintf('\n%s (Aggregated Ensemble) => CV1 [%g, %g, regression]:\t', ...
                                                    strout,i,j); 
                                        end
                                    end
                                    % Compute ensemble diversity criterion
                                    % only if more than 1 base learner is
                                    % present
                                    OUT.TrDiv(i,j,curclass) = feval( DivFunc, OUT.TrHT{i,j,curclass}, TrL{curclass} );
                                    OUT.CVDiv(i,j,curclass) = feval( DivFunc, OUT.CVHT{i,j,curclass}, CVL{curclass} );
                                     % Evaluate ensemble performance
                                    OUT.TrHDperf(i,j,curclass) = nk_EnsPerf( OUT.TrHD{i,j,curclass}, TrL{curclass} );
                                    OUT.TrHTperf(i,j,curclass) = nk_EnsPerf( OUT.TrHT{i,j,curclass}, TrL{curclass} );       
                                    OUT.CVHDperf(i,j,curclass) = nk_EnsPerf( OUT.CVHD{i,j,curclass}, CVL{curclass} );
                                    OUT.CVHTperf(i,j,curclass) = nk_EnsPerf( OUT.CVHT{i,j,curclass}, CVL{curclass} );
                                else
                                    OUT.TrHDperf(i,j,curclass) = OUT.tr{i,j,curclass}(kInd);
                                    OUT.TrHTperf(i,j,curclass) = OUT.tr{i,j,curclass}(kInd);
                                    OUT.CVHDperf(i,j,curclass) = OUT.ts{i,j,curclass}(kInd);
                                    OUT.CVHTperf(i,j,curclass) = OUT.ts{i,j,curclass}(kInd);
                                end                             

                                tkIndCat = [tkIndCat kInd'];
                                if ns > 1
                                    tSubSpaces{i,j,curclass}    = IN.F{i,j,curclass}(:,OUT.kxVec{i,j,curclass}(kIndN{curclass}));
                                    tF{i,j,curclass}            = kInd; 
                                    tWeights{i,j,curclass}      = W{curclass}(kIndN{curclass});
                                    tModels{i,j,curclass}       = OUT.mdl{i,j,curclass}(kIndN{curclass});

                                    if W2AVAIL
                                        tW2{i,j,curclass}       = OUT.w2{i,j,curclass}(kIndN{curclass});
                                        tMd{i,j,curclass}       = OUT.Md{i,j,curclass}(kIndN{curclass});
                                        tMm{i,j,curclass}       = OUT.Mm{i,j,curclass}(kIndN{curclass});
                                    end
                                else
                                    tSubSpaces{i,j,curclass}    = IN.F{i,j}(:,OUT.kxVec{i,j,curclass}(kIndN{curclass}));
                                    tF{i,j,curclass}            = kInd; 
                                    tWeights{i,j}               = W{curclass}(kIndN{curclass});
                                    tModels{i,j,curclass}       = OUT.mdl{i,j,curclass}(kIndN{curclass});

                                    if W2AVAIL
                                        tW2{i,j,curclass}       = OUT.w2{i,j,curclass}(kIndN{curclass});
                                        tMd{i,j,curclass}       = OUT.Md{i,j,curclass}(kIndN{curclass});
                                        tMm{i,j,curclass}       = OUT.Mm{i,j,curclass}(kIndN{curclass});
                                    end
                                end

                        otherwise
                                %% OPTIMIZED ENSEMBLE CONSTRUCTION 
                                 % (Recursive classifier elimination of forward classifier construction)                        
                                if ~MULTI.flag || ~MULTI.train % Perform binary optimization

                                    switch Param.EnsembleStrategy.DataType

                                        case 1 % Use CV1 training data for ensemble construction
                                            % Only for binary optimization: Use CV1-training data to optimize ensemble components
                                            Lx = TrL{curclass};
                                            
                                            switch Metric
                                                case 1 % Targets
                                                    
                                                    Px = OUT.TrHT{i,j,curclass};
                                                    [tkInd, Hx_perf, Hx, Hx_div ] = nk_BuildEnsemble(Px, Lx, Param.EnsembleStrategy,[], ngroups);
                                                    OUT.TrHTperf(i,j,curclass) = Hx_perf; OUT.TrHT{i,j,curclass} = Hx;
                                                    OUT.TrHD{i,j,curclass} = OUT.TrHD{i,j,curclass}(:,tkInd);
                                                    
                                                case 2 % Decision values
                                                    
                                                    Px = OUT.TrHD{i,j,curclass};  
                                                    [tkInd, Hx_perf, Hx, Hx_div ] = nk_BuildEnsemble(Px, Lx, Param.EnsembleStrategy,[], ngroups);
                                                    OUT.TrHDperf(i,j,curclass) = Hx_perf; OUT.TrHD{i,j,curclass} = Hx;
                                                    OUT.TrHT{i,j,curclass} = OUT.TrHT{i,j,curclass}(:,tkInd);
                                            end
                                            
                                            OUT.TrDiv(i,j,curclass)     = Hx_div;
                                            % Extract base classifiers from
                                            % CV1 test ensemble according to tkInd 
                                            OUT.CVHD{i,j,curclass}       = OUT.CVHD{i,j,curclass}(:,tkInd); 
                                            OUT.CVHT{i,j,curclass}       = OUT.CVHT{i,j,curclass}(:,tkInd);
                                            % Update performance of optimized CV1 test data ensemble
                                            OUT.CVHDperf(i,j,curclass)   = nk_EnsPerf( OUT.CVHD{i,j,curclass}, CVL{curclass} );
                                            OUT.CVHTperf(i,j,curclass)   = nk_EnsPerf( OUT.CVHT{i,j,curclass}, CVL{curclass} );
                                            % Compute entropy of optimized CV1 test data ensemble
                                            OUT.CVDiv(i,j,curclass) = feval( DivFunc ,OUT.CVHT{i,j,curclass}, CVL{curclass} );
                                            
                                        case 2 % Use CV1 test data for ensemble construction
                                            % Only for binary optimization: Use CV1-test data to optimize ensemble components
                                            Lx = CVL{curclass};
                                            
                                            switch Metric
                                                case 1
                                            
                                                    Px = OUT.CVHT{i,j,curclass}; 
                                                    [tkInd, Hx_perf, Hx, Hx_div ] = nk_BuildEnsemble(Px, Lx, Param.EnsembleStrategy,[], ngroups);
                                                    OUT.CVHTperf(i,j,curclass) = Hx_perf; OUT.CVHT{i,j,curclass} = Hx;
                                                    OUT.CVHD{i,j,curclass} = OUT.CVHD{i,j,curclass}(:,tkInd);
                                                    
                                                case 2
                                                    
                                                    Px = OUT.CVHD{i,j,curclass}; 
                                                    [ tkInd, Hx_perf, Hx, Hx_div ] = nk_BuildEnsemble(Px, Lx, Param.EnsembleStrategy,[],ngroups);
                                                    OUT.CVHDperf(i,j,curclass)  = Hx_perf; OUT.CVHD{i,j,curclass}      = Hx;  
                                                    OUT.CVHT{i,j,curclass} = OUT.CVHT{i,j,curclass}(:,tkInd);
                                            end
                                            
                                            OUT.CVDiv(i,j,curclass)     = Hx_div;    
                                            % Extract base classifiers from
                                            % CV1 training ensemble according to tkInd 
                                            OUT.TrHD{i,j,curclass}       = OUT.TrHD{i,j,curclass}(:,tkInd); 
                                            OUT.TrHT{i,j,curclass}       = OUT.TrHT{i,j,curclass}(:,tkInd);
                                            % Update performance of optimized CV1 training data ensemble
                                            OUT.TrHDperf(i,j,curclass)   = nk_EnsPerf( OUT.TrHD{i,j,curclass}, TrL{curclass} );
                                            OUT.TrHTperf(i,j,curclass)   = nk_EnsPerf( OUT.TrHT{i,j,curclass}, TrL{curclass} );
                                            % Compute entropy of optimized CV1 training data ensemble
                                            OUT.TrDiv(i,j,curclass) = feval( DivFunc, OUT.TrHT{i,j,curclass}, TrL{curclass});
                                            
                                        case 3 % Use CV1 training & test data for ensemble construction
                                             
                                            Lx = [TrL{curclass}; CVL{curclass}];
                                            %Ix = [ones(size(TrL{curclass},1),1); 2*ones(size(CVL{curclass},1),1)];
                                            
                                            switch Metric
                                                case 1 
                                                    Px = [ OUT.TrHT{i,j,curclass}; OUT.CVHT{i,j,curclass} ]; 
                                                    tkInd = nk_BuildEnsemble(Px, Lx, Param.EnsembleStrategy, [], ngroups);
                                                    
                                                case 2
                                                    Px = [ OUT.TrHD{i,j,curclass}; OUT.CVHD{i,j,curclass} ]; 
                                                    tkInd = nk_BuildEnsemble(Px, Lx, Param.EnsembleStrategy,[],ngroups);
                                            end
                                           
                                            OUT.TrHT{i,j,curclass}      = OUT.TrHT{i,j,curclass}(:,tkInd);
                                            OUT.CVHT{i,j,curclass}      = OUT.CVHT{i,j,curclass}(:,tkInd);
                                            OUT.TrHD{i,j,curclass}      = OUT.TrHD{i,j,curclass}(:,tkInd);
                                            OUT.CVHD{i,j,curclass}      = OUT.CVHD{i,j,curclass}(:,tkInd);
                                            OUT.TrDiv(i,j,curclass)     = feval( DivFunc, OUT.TrHT{i,j,curclass}, TrL{curclass});
                                            OUT.CVDiv(i,j,curclass)     = feval( DivFunc, OUT.CVHT{i,j,curclass}, CVL{curclass});
                                            OUT.TrHDperf(i,j,curclass)  = nk_EnsPerf( OUT.TrHD{i,j,curclass}, TrL{curclass} );
                                            OUT.TrHTperf(i,j,curclass)  = nk_EnsPerf( OUT.TrHT{i,j,curclass}, TrL{curclass} );
                                            OUT.CVHDperf(i,j,curclass)  = nk_EnsPerf( OUT.CVHD{i,j,curclass}, CVL{curclass} );
                                            OUT.CVHTperf(i,j,curclass)  = nk_EnsPerf( OUT.CVHT{i,j,curclass}, CVL{curclass} );
                                            
                                    end
                                    % Is this correct?
                                    tInd = false(1,length(kIndN{curclass})); tInd(tkInd)=true;
                                    tkIndCat = [tkIndCat tInd];
                                    % Update feature subspace mask indices
                                    kIndN{curclass} = kIndN{curclass}(tkInd);
                                    % Update feature subspaces
                                    tSubSpaces{i,j,curclass} = IN.F{i,j,curclass}(:,OUT.kxVec{i,j,curclass}(kIndN{curclass}));
                                    % Update feature subspace masks
                                    tF{i,j,curclass} = false(size(OUT.F{i,j,curclass}));
                                    tF{i,j,curclass}(kIndN{curclass}) = true;
                                    % Update feature subspace weights
                                    tWeights{i,j,curclass} = W{curclass}(tkInd);
                                    tModels{i,j,curclass}  = OUT.mdl{i,j,curclass}(kIndN{curclass});

                                    if W2AVAIL
                                        tW2{i,j,curclass}       = OUT.w2{i,j,curclass}(kIndN{curclass});
                                        tMd{i,j,curclass}       = OUT.Md{i,j,curclass}(kIndN{curclass});
                                        tMm{i,j,curclass}       = OUT.Mm{i,j,curclass}(kIndN{curclass});
                                    end
                                end     
                    end
                end

                %% MULTI-GROUP OPTIMIZATION                              
                if MULTI.flag

                    % Hardcode if ECOC will be used:
                    if MULTI.method == 2, mTr = sign(mTr); mCV = sign(mCV); end

                    switch EnsType

                        case 0
                            %% MULTI-GROUP AGGREGATED ENSEMBLE (no classifier selection)
                            % Compute CV1-training and test data
                            % performance using simply the aggregated
                            % multi-group ensembles:
                            [OUT.mTrPerf(i,j), OUT.mTrPred{i,j}] = nk_MultiEnsPerf(mTr, sign(mTr), IN.Y.mTrL{i,j}(:,MULTILABEL.curdim), mC, ngroups);
                            [OUT.mCVPerf(i,j), OUT.mCVPred{i,j}] = nk_MultiEnsPerf(mCV, sign(mCV), IN.Y.mCVL{i,j}(:,MULTILABEL.curdim), mC, ngroups);
                            % Compute CV1-training and test data
                            % ensemble entropy as the mean of
                            % dichotomizers' entropies
                            tTrDiv = zeros(IN.nclass,1); tCVDiv = tTrDiv;
                            if ~isempty(DivFunc)
                                for qcurclass = 1 : IN.nclass
                                    indC = mC==qcurclass;
                                    tTrDiv(curclass) = feval( DivFunc, mTr(:,indC), TrL{qcurclass} );
                                    tCVDiv(curclass) = feval( DivFunc, mCV(:,indC), CVL{qcurclass} );
                                end
                                OUT.mTrDiv(i,j)     = sum(tTrDiv) / IN.nclass;
                                OUT.mCVDiv(i,j)     = sum(tCVDiv) / IN.nclass;
                            end
                            % Subspace & models don't have to be
                            % updated, this has already been done in
                            % the binary classifier section

                        otherwise
                            %% OPTIMIZED MULTI-GROUP ENSEMBLE CONSTRUCTION 
                            if MULTI.train

                                % Perform multi-class optimization
                                switch Param.EnsembleStrategy.DataType
                                    case 1 % CV1-Training Data
                                        [tkInd, ...
                                            OUT.mTrPerf(i,j), ...
                                            dum, ...
                                            OUT.mTrDiv(i,j), ...
                                            tkIndCat, ...
                                            OUT.mTrPred{i,j}] = ...
                                                        nk_BuildEnsemble(mTr, IN.Y.mTrL{i,j}(:,MULTILABEL.curdim), Param.EnsembleStrategy, mC, ngroups);

                                        % Compute CV1-test data multi-group performance:
                                        [OUT.mCVPerf(i,j), OUT.mCVPred{i,j}] = nk_MultiEnsPerf(mCV(:,tkIndCat), sign(mCV(:,tkIndCat)), IN.Y.mCVL{i,j}(:,MULTILABEL.curdim), mC(tkIndCat), ngroups);
                                        % Compute CV1-test data ensemble entropy:
                                        tDiv = zeros(IN.nclass,1);
                                        for qcurclass = 1 : IN.nclass, 
                                            tDiv(curclass) = feval( DivFunc, mCV(:,tkInd{qcurclass}), CVL{qcurclass});
                                        end
                                        OUT.mCVDiv(i,j) = sum(tDiv) / IN.nclass; 
                                    case 2 % CV1-Test Data
                                        [tkInd, ...
                                            OUT.mCVPerf(i,j), ...
                                            dum, ...
                                            OUT.mCVDiv(i,j), ...
                                            tkIndCat,...
                                            OUT.mCVPred{i,j}] = ...
                                                        nk_BuildEnsemble(mCV, IN.Y.mCVL{i,j}(:,MULTILABEL.curdim), Param.EnsembleStrategy, mC, ngroups);
                                        % Compute CV1-training data multi-group performance:
                                        [OUT.mTrPerf(i,j), OUT.mTrPred{i,j}] = ...
                                            nk_MultiEnsPerf(mTr(:,tkIndCat), sign(mTr(:,tkIndCat)), IN.Y.mTrL{i,j}(:,MULTILABEL.curdim), mC(tkIndCat), ngroups);
                                        % Compute CV1-training data ensemble entropy:
                                        tDiv = zeros(IN.nclass,1);
                                        for qcurclass = 1 : IN.nclass, 
                                            tDiv(curclass) = feval( DivFunc, mTr(:,tkInd{qcurclass}), TrL{qcurclass} );
                                        end
                                        OUT.mTrDiv(i,j) = sum(tDiv) / IN.nclass; 
                                end

                                % Apply tkInd to dichotomization info
                                for curclass = 1 : IN.nclass
                                    % Update feature mask indices with
                                    % optimization info:
                                    kIndN{curclass}               = kIndN{curclass}(tkInd{curclass});
                                    % Apply this to the training data:
                                    OUT.TrHD{i,j,curclass}          = OUT.TrHD{i,j,curclass}(:,tkInd{curclass}); 
                                    OUT.TrHT{i,j,curclass}          = OUT.TrHT{i,j,curclass}(:,tkInd{curclass});
                                    OUT.TrHDperf(i,j,curclass)      = nk_EnsPerf( OUT.TrHD{i,j,curclass}, TrL{curclass} );
                                    OUT.TrHTperf(i,j,curclass)      = nk_EnsPerf( OUT.TrHT{i,j,curclass}, TrL{curclass} );
                                    OUT.TrDiv(i,j,curclass)         = feval( DivFunc, OUT.TrHT{i,j,curclass}, TrL{curclass} );
                                    % ... to the CV1-test data:
                                    OUT.CVHD{i,j,curclass}          = OUT.CVHD{i,j,curclass}(:,tkInd{curclass} ); 
                                    OUT.CVHT{i,j,curclass}          = OUT.CVHT{i,j,curclass}(:,tkInd{curclass} );
                                    OUT.CVHDperf(i,j,curclass)      = nk_EnsPerf( OUT.CVHD{i,j,curclass}, CVL{curclass} );
                                    OUT.CVHTperf(i,j,curclass)      = nk_EnsPerf( OUT.CVHT{i,j,curclass}, CVL{curclass} );
                                    OUT.CVDiv(i,j,curclass)         = feval( DivFunc, OUT.CVHT{i,j,curclass}, CVL{curclass} );
                                    % ... and to the feature subspace info:
                                    if ns > 1
                                        tSubSpaces{i,j,curclass}    = IN.F{i,j,curclass}(:,OUT.kxVec{i,j,curclass}(kIndN{curclass}));
                                    else
                                        tSubSpaces{i,j,curclass}    = IN.F{i,j}(:,OUT.kxVec{i,j,curclass}(kIndN{curclass}));
                                    end
                                    if nc > 1
                                        tF{i,j,curclass}            = false(size(OUT.F{i,j,curclass}));
                                    else
                                        tF{i,j,curclass}            = false(size(OUT.F{i,j}));
                                    end
                                    tF{i,j,curclass}(kIndN{curclass}) = true;
                                    tWeights{i,j,curclass}          = W{curclass}(tkInd{curclass});
                                    tModels{i,j,curclass}           = OUT.mdl{i,j,curclass}(kIndN{curclass});
                                    if W2AVAIL
                                        tW2{i,j,curclass}       = OUT.w2{i,j,curclass}(kIndN{curclass});
                                        tMd{i,j,curclass}       = OUT.Md{i,j,curclass}(kIndN{curclass});
                                        tMm{i,j,curclass}       = OUT.Mm{i,j,curclass}(kIndN{curclass});
                                    end
                                end
                            else
                                % Simpy extract the multi-group
                                % ensemble by applying tkIndCat to the
                                % concated dichotomization data
                                tkIndCat = logical(tkIndCat);
                                mTr = mTr(:,tkIndCat); mCV = mCV(:,tkIndCat); mC = mC(tkIndCat);
                                % CV1-training data ensemble
                                % performance and entropy
                                [OUT.mTrPerf(i,j), OUT.mTrPred{i,j}] = nk_MultiEnsPerf( mTr, mTr, IN.Y.mTrL{i,j}, mC , ngroups);
                                OUT.mTrDiv(i,j) = feval( DivFunc, mTr, TrL{curclass} );
                                % ... the same with CV1-test data:
                                [OUT.mCVPerf(i,j), OUT.mCVPred{i,j}] = nk_MultiEnsPerf( mCV, mCV, IN.Y.mCVL{i,j}, mC, ngroups );
                                OUT.mCVDiv(i,j) = feval( DivFunc, mCV, CVL{curclass} );
                            end

                    end
                end
        end
    end
end

if EnsType ~=9
    % Update IN and OUT
    IN.F        = tSubSpaces;
    OUT.F       = tF;
    OUT.Weights = tWeights;
    OUT.mdl     = tModels;
    if W2AVAIL
        OUT.w2  = tW2;
        OUT.Md  = tMd;
        OUT.Mm  = tMm;
    end
end

end
