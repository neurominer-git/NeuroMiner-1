%==========================================================================
% [IN, OUT] = FoldPerm(IN, OUT, strout, fRFE, fFull, fKX, LoopParam)
%==========================================================================
% MAIN NM TRAINING MODULE 
% FoldPerm manages the logic of model training and optimization at the
% inner cross-validation cycle. It main operational mode is to generate 
% the models using the CV1 training data and to apply them to CV1 test data 
% (with or without filtering or wrapping). Alternatively, it can bypass 
% this process and train models using the entire CV1 data partition in case
% no hyperparameter optimization is needed or in case models need to be re-
% trained after model optimization
%
% Inputs:
% -----------
% IN        :   Parameter structure governing the model optimization process
% OUT       :   Results obtained from the model optimization process
% strout    :   String describing the operational mode of FoldPerm
% fRFE      :   flag indicating whether FoldPerm needs to conduct
%               wrapper-based feature selection
% fFull     :   flag indicating whether FoldPerm is using the entire CV1
%               partition data
% fKX       :   flag/integer telling FoldPerm how to conduct subspace-based
%               filtering
%
% Outputs: IN & OUT, see above
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%(c) Nikolaos Koutsouleris, 03/2020

function [IN, OUT] = FoldPerm(IN, OUT, strout, fRFE, fFull, RetrainImmediate, fKX, LoopParam)

global VERBOSE CV MODEFL RAND MULTILABEL MULTI RFE CVPOS

RF = []; fMULTI = false; VI = []; 
if MULTI.flag ==1 && ... 
        MULTI.train == 1 && ...
        ( isfield(RFE.Wrapper.GreedySearch,'MultiClassOptimization') && ...
            RFE.Wrapper.GreedySearch.MultiClassOptimization)
    fMULTI = true; 
end

%% Setup cell array containers
if ~exist('OUT', 'var') || isempty(OUT)
    
    OUT.tr      = cell( IN.nperms, IN.nfolds, IN.nclass ); % Training performance
    OUT.mdl     = cell( IN.nperms, IN.nfolds, IN.nclass ); % Model structures
    OUT.Trtargs = cell( IN.nperms, IN.nfolds, IN.nclass ); % Predicted target labels for CV1 training data
    OUT.Trdecs  = cell( IN.nperms, IN.nfolds, IN.nclass ); % Decision values / probabilities for CV1 training data
    OUT.kxVec   = cell( IN.nperms, IN.nfolds, IN.nclass ); % Subspace Stepping
    OUT.ts      = cell( IN.nperms, IN.nfolds, IN.nclass ); % Test performance
    OUT.rf      = cell( IN.nperms, IN.nfolds, IN.nclass ); % Recursive feature elimination structures (unused)
    OUT.w2      = cell( IN.nperms, IN.nfolds, IN.nclass ); % |w| for linear SVM kernels
    OUT.Md      = cell( IN.nperms, IN.nfolds, IN.nclass ); % distance to hyperplane for linear SVM kernels
    OUT.Mm      = cell( IN.nperms, IN.nfolds, IN.nclass ); % normalized margin for linear SVM kernels
    OUT.CVtargs = cell( IN.nperms, IN.nfolds, IN.nclass ); % Predicted target labels for CV1 test data
    OUT.CVdecs  = cell( IN.nperms, IN.nfolds, IN.nclass ); % Decision values / probabilities for CV1 test data
    OUT.featnum = cell( IN.nperms, IN.nfolds, IN.nclass, IN.nvar ); % Number of feature
    
end
nc = size(IN.F,3); % number of binary comparisons in IN.F (features)

if exist('LoopParam','var') && ~isempty(LoopParam)
    PermVec = LoopParam.PermVec;
    FoldVec = LoopParam.FoldVec;
    ClassVec = LoopParam.ClassVec;
else
    PermVec = 1:IN.nperms;
    FoldVec = 1:IN.nfolds; 
    ClassVec = 1:IN.nclass;
end

PermNum = numel(PermVec);
FoldNum = numel(FoldVec);
ClassNum = numel(ClassVec);
kx = 1;
CVPOS.fFull = fFull;
for ii=1:PermNum % Loop through CV1 permutations

    for jj=1:FoldNum % Loop through CV1 folds
        
        % Initialize variables for partition
        i           = PermVec(ii); 
        j           = FoldVec(jj);
        CVPOS.CV1p = i;
        CVPOS.CV1f = j;
        modelTrL    = cell(1,ClassNum);
        Ymodel      = cell(1,ClassNum);
        Ytrain      = cell(1,ClassNum);
        Ytest       = cell(1,ClassNum);
        tTrL        = cell(1,ClassNum);
        tCVL        = cell(1,ClassNum);
        kFea        = zeros(1,ClassNum);
        lFea        = zeros(1,ClassNum);
        Fx          = cell(1,ClassNum);
        cvts_fl     = false(1,ClassNum);
        k           = zeros(1,ClassNum);
        Fk          = cell(ClassNum, IN.nvar); 
        
%       PART I
%       Loop through dichotomizers and prepare feature spaces and labels for training
%       -----------------------------------------------------------------------------
        for ccurclass=1:ClassNum 
            
            % Initialize variable for dichotomizer
            curclass = ClassVec(ccurclass);
            tFea = zeros(curclass, IN.nvar); 
            tCV = cell(1,IN.nvar); 
            tTr = cell(1,IN.nvar);
            modelTr = cell(1,IN.nvar);
            CVInd = IN.Y.CVInd{i,j}{curclass};
            TrInd = IN.Y.TrInd{i,j}{curclass};
            
            %% Loop through variates and extract multi-modal training data cell array
            for v=1:IN.nvar 
                
                % Get feature subspace size
                if nc > 1 % more than one binary classifier
                    try 
                        tFea(curclass,v) = size(IN.F{i,j,curclass,v},2); 
                    catch
                        tFea(curclass,v) =  size(IN.F{i,j,curclass,v},1); 
                    end
                    Fk{curclass,v} = IN.F{i,j,curclass,v}; 
                else
                    try % a binary classifier (a multi-group classifier, not impl.)
                        tFea(curclass,v) = size(IN.F{i,j,v},2); 
                    catch
                        tFea(curclass,v) =  size(IN.F{i,j,v},1); 
                    end
                    Fk{curclass,v} = IN.F{i,j,v}; 
                end
                
                % Determine training data & validation data
                if ~iscell(IN.Y.Tr{i,j,1})
                    modelTr{v} = IN.Y.Tr{i,j,v}(TrInd,:);
                    tTr{v} = IN.Y.Tr{i,j,v};
                    if fFull, modelTr{v} = [modelTr{v}; IN.Y.CV{i,j,v}(CVInd,:) ]; end
                    tCV{v} = IN.Y.CV{i,j,v};
                else
                    modelTr{v} = IN.Y.Tr{i,j,v}{curclass}(TrInd,:);
                    tTr{v} = IN.Y.Tr{i,j,v}{curclass};
                    if fFull, modelTr{v} = [modelTr{v}; IN.Y.CV{i,j,v}{curclass}(CVInd,:) ]; end
                    tCV{v} = IN.Y.CV{i,j,v}{curclass};
                end
              
            end
            
            %% Get samples sizes
            kSubjTr = size(tTr{1},1); 
            kSubjCV = size(tCV{1},1);
            
            %% Determine labels for learning process
            modelTrL{curclass} = IN.Y.TrL{i,j}{curclass}(:,MULTILABEL.curdim);
            tTrL{curclass} = zeros(size(tTr{1},1),1); 
            tTrL{curclass}(TrInd) = IN.Y.TrL{i,j}{curclass}(:,MULTILABEL.curdim);
            tCVL{curclass} = zeros(size(tCV{1},1),1); 
            tCVL{curclass}(CVInd) = IN.Y.CVL{i,j}{curclass}(:,MULTILABEL.curdim);
            if fFull, 
                modelTrL{curclass} = [modelTrL{curclass}; IN.Y.CVL{i,j}{curclass}(:,MULTILABEL.curdim)]; 
            end
            
            %% Define feature subspace counter ...
            % according to maximum
            % feature subspace size across variates. Compute subspace
            % stepping for loop.
            % Stepping for feature subspace search
            lFea(curclass) = max(tFea(curclass,:)); 
            
            if fKX, kx = ceil((lFea(curclass) / 100) * fKX); end
            
            OUT.kxVec{i,j,curclass} = kx:kx:lFea(curclass);
            if isempty(OUT.kxVec{i,j,curclass})
                OUT.kxVec{i,j,curclass} = 1;
            else
                if OUT.kxVec{i,j,curclass}(end) < lFea(curclass), 
                    OUT.kxVec{i,j,curclass} = [ OUT.kxVec{i,j,curclass} lFea(curclass) ]; 
                end; 
            end
            kFea(curclass) = length(OUT.kxVec{i,j,curclass});
            cvts_fl(curclass) = kFea(curclass) ~= size(OUT.ts{i,j,curclass},1);
            
            %% Initialize arrays for current multi-dimensional cell pointer
            OUT.tr{i,j,curclass}         = zeros( kFea(curclass), 1 );
            if isempty(OUT.mdl{i,j,curclass})
                OUT.mdl{i,j,curclass}        = cell( kFea(curclass), 1 );
            end
            OUT.Trtargs{i,j,curclass}    = zeros( kSubjTr, kFea(curclass) );
            OUT.Trdecs{i,j,curclass}     = zeros( kSubjTr, kFea(curclass) );
            
            if ~fFull || RetrainImmediate || any(cvts_fl)
                OUT.ts{i,j,curclass}     = zeros( kFea(curclass),1 );
                OUT.CVtargs{i,j,curclass}= zeros( kSubjCV, kFea(curclass) );
                OUT.CVdecs{i,j,curclass} = zeros( kSubjCV, kFea(curclass) );
                OUT.rf{i,j,curclass}     = cell( kFea(curclass), 1 );
                OUT.w2{i,j,curclass}     = zeros( kFea(curclass), 1 );
                OUT.Md{i,j,curclass}     = zeros( kFea(curclass), 1 );
                OUT.Mm{i,j,curclass}     = zeros( kFea(curclass), 1 );
            end
            
            for kT=1:kFea(curclass)
                
                k(curclass) = OUT.kxVec{i,j,curclass}(kT);
                
                %% Extract training and test data according to current feature subspace mask Fk
                [Ymodel{curclass}{k(curclass) }, Fx{curclass}{k(curclass)}] = nk_ExtractFeatures(modelTr, Fk(curclass,:), [], k(curclass) );
                Ytrain{curclass}{k(curclass) } = nk_ExtractFeatures(tTr, Fk(curclass,:), [], k(curclass) );
                Ytest{curclass}{k(curclass) }  = nk_ExtractFeatures(tCV, Fk(curclass,:), [], k(curclass) );
            end
        end
        
%       PART II
%       Loop through dichotomizers and train models using prepared feature spaces
%       Eventually use multi-group optimization! 
%       --------------------------------------------------------------------------
        if fMULTI && fRFE
         
            tMTrL   = IN.Y.mTrL{i,j};
            tMCVL   = IN.Y.mCVL{i,j};
            indkX = 1;

            % Assign optimization results to dichotomizers                
            for kT=1:kFea(curclass)
                
                tYmodel = cell(1, ClassNum);
                tYtrain = cell(1, ClassNum);
                tYtest   = cell(1, ClassNum);
                %tModelL = cell(1, Classnum);
                %tTrL    = cell(1, ClassNum);
                
                 % use wrapper-based optimization (binary models / regression)
                if IN.nvar < 2 % Univariate case
                    OUT.featout{i,j,curclass} = zeros(size(IN.F{i,j,curclass},1),kFea(curclass)); 
                else % Multiple variate case
                    OUT.featout = cell(IN.nvar,1);
                    for v=1:IN.nvar 
                        OUT.featout{v}{i,j,curclass} = zeros(size(IN.F{i,j,curclass,v},1),tFea(curclass,v)); 
                    end
                end
                
                for ccurclass=1:ClassNum 

                    curclass            = ClassVec(ccurclass);
                    tYmodel{curclass}   = Ymodel{curclass}{k(curclass)};
                    tYtrain{curclass}   = Ytrain{curclass}{k(curclass)};
                    tYtest{curclass}    = Ytest{curclass}{k(curclass)};
                    
                end
                % Perform multi-class wrapper-based optimization
                [RF, model]  = nk_MLOptimizer_WrapperMulti(tYmodel, tYtrain, modelTrL, tTrL, tMTrL, tYtest, tCVL, tMCVL, IN.ngroups, IN.Ps, []);
                
                 % Assign optimization results to dichotomizers
                for ccurclass=1:ClassNum 
                    curclass = ClassVec(ccurclass);
                    
                    % Update data and feature masks according to wrapper-based                                   
                    [OUT, IN, Ytrain, Ytest, Ymodel] = assignRF2data(OUT, RF, IN, Ytrain, Ytest, Ymodel, Fx, i, j, curclass, k, indkX, v);

                    % Compute model performance in training and test data
                    OUT = getperf2out(OUT, model, Ymodel, Ytrain, tTrL, Ytest, tCVL, fFull, RetrainImmediate, cvts_fl, i, j, curclass, indkX, k);

                    % Store number of feeatures in OUT if no wrapper-based
                    % mask available or detected
                    if isempty(RF) || ~RF.found
                        for v=1:IN.nvar 
                             OUT.featnum{i,j,curclass,v}(indkX) = sum(Fk{curclass, v}(:,indkX)~=0);
                        end
                    end

                    % Assign training and test results to OUT structure
                    OUT = assign2out( OUT, RetrainImmediate, fFull, cvts_fl, i , j, curclass);
                    
                end
                indkX = indkX+1;    
             
            end   

        else
            
            % Optimize dichtomizers separately
            for ccurclass=1:ClassNum 

                curclass = ClassVec(ccurclass);

                % Loop through every (or every kx, to save run time) feature subspace
                indkX = 1;

                for kT=1:kFea(curclass)

                    nfeats = sum(any(Fk{1}(:,k),2));
                    if VERBOSE, 
                        switch MODEFL
                            case 'classification'
                               if RAND.Decompose ~= 9
                                    fprintf('\n%s => CV1 [%g, %g, %s, Suspace %g/%g (%g feats)]:', strout, ...
                                        i, j, CV.class{1,1}{curclass}.groupdesc, k(curclass) , lFea(curclass), nfeats)
                               else
                                   fprintf('\n%s => CV1 [%g, %g, Multi-Group, Subspace %g/%g (%g feats)]:', strout, ...
                                        i, j, k(curclass) , lFea(curclass), nfeats)
                               end
                            case 'regression'
                               fprintf('\n%s => CV1 [%g, %g, Regression, Subspace %g/%g (%g feats)]:', strout, ...
                                        i, j, k(curclass) , lFea(curclass), nfeats)
                        end
                    end
                    if isfield(IN.Y,'VI'), VI =  IN.Y.VI{ i,j}{curclass}; end
                    
                    %% Train model(s)
                    if ~fRFE

                       % Train algorithm without wrapper
                        [~, model]= nk_GetParam2(Ymodel{curclass}{k(curclass)}, modelTrL{curclass}, IN.Ps{curclass}, 1, VI);
                        
                    else
                        
                        % use wrapper-based optimization (binary models / regression)
                        if IN.nvar < 2 % Univariate case
                            OUT.featout{i,j,curclass} = zeros(size(IN.F{i,j,curclass},1),kFea(curclass)); 
                        else % Multiple variate case
                            OUT.featout = cell(IN.nvar,1);
                            for v=1:IN.nvar 
                                OUT.featout{v}{i,j,curclass} = zeros(size(IN.F{i,j,curclass,v},1),tFea(curclass,v)); 
                            end
                        end

                        % NOTE: RFE currently supports only univariate
                        % data!!!
                        [RF, model]  = nk_MLOptimizer_Wrapper(Ymodel{curclass}{k(curclass)}, modelTrL{curclass}, ...
                                                              Ytest{curclass}{k(curclass)}, tCVL{curclass}, IN.Ps{curclass}, []);

                        % Update data and feature masks according to wrapper-based                                   
                        [OUT, IN, Ytrain, Ytest, Ymodel] = assignRF2data(OUT, RF, IN, Ytrain, Ytest, Ymodel, Fx, i, j, curclass, k, indkX, v);
                           
                    end
                    % Compute model performance in training and test data
                    OUT = getperf2out(OUT, model, Ymodel, Ytrain, tTrL, Ytest, tCVL, fFull, RetrainImmediate, cvts_fl, i, j, curclass, indkX, k);
                    
                    % Store number of feeatures in OUT if no wrapper-based
                    % mask available or detected
                    if isempty(RF) || ~RF.found
                        for v=1:IN.nvar 
                             OUT.featnum{i,j,curclass,v}(indkX) = sum(Fk{curclass, v}(:,indkX)~=0);
                        end
                    end
                    % Assign training and test results to OUT structure
                    OUT = assign2out( OUT, RetrainImmediate, fFull, cvts_fl, i , j, curclass);
                    
                    indkX = indkX+1;
                end
            end
        end
    end
end

% _________________________________________________________________________
function [OUT, IN, Ytrain, Ytest, Ymodel] = assignRF2data(OUT, RF, IN, ...
                                                    Ytrain, Ytest, Ymodel, ...
                                                    Fx, i, j, curclass, ...
                                                    k, indkX, v)

if RF.found
    % Transfer new subspace features to SubSets 
    indFx                                   = find(Fx{curclass}{k(curclass)}); 
    if iscell(RF.FeatureIndex)
        FI = RF.FeatureIndex{curclass};
        if isempty(FI), error('Wrapper of model #%g did not returned any features! Relax your wrapper settings', curclass); end
    else
        FI = RF.FeatureIndex;
        if isempty(FI), error('Wrapper did not returned any features! Relax your wrapper settings'); end
    end
    indFx                                   = indFx(FI);
    indFx0                                  = zeros( size( Fx{curclass}{k(curclass)}));  
    indFx0(indFx)                           = IN.F{i,j,curclass,v}(indFx,k(curclass) );
    % Extract RFE subspace from Ytrain and Ytest
    Ytrain{curclass}{k(curclass)}           = Ytrain{curclass}{k(curclass)}(:,FI);
    Ytest{curclass}{k(curclass)}            = Ytest{curclass}{k(curclass)}(:,FI);
    Ymodel{curclass}{k(curclass)}           = Ymodel{curclass}{k(curclass)}(:,FI);
    IN.F{i,j,curclass,v}(:,k(curclass) )    = indFx0;
    OUT.featout{i,j,curclass}(:,indkX)      = indFx0;
    OUT.featnum{i,j,curclass}(indkX)        = sum(OUT.featout{i,j,curclass}(:,indkX)~=0); 
    
end

% _________________________________________________________________________
function OUT = getperf2out(OUT, model, Ymodel, Ytrain, tTrL, Ytest, tCVL, ...
                            fFull, RetrainImmediate, cvts_fl, ...
                            i, j, curclass, indkX, k)
global VERBOSE

if iscell(model)
    iModel = model{curclass};
else
    iModel = model;
end

 %% Apply trained algorithm to CV1 test data
[OUT.tr{i,j,curclass}(indkX), ...
    OUT.Trtargs{i,j,curclass}(:,indkX) , ...
    OUT.Trdecs{i,j,curclass}(:,indkX), iModel ] = nk_GetTestPerf(Ytrain{curclass}{k(curclass)}, tTrL{curclass}, ...
                                                         [], iModel, Ymodel{curclass}{k(curclass)});
if isnan(OUT.tr{i,j,curclass}(indkX))
    warning('Non-finite performance measures found in CV1 training data')
end
if VERBOSE, fprintf('\tTr = %1.2f', OUT.tr{i,j,curclass}(indkX)); end
if ~fFull || RetrainImmediate || cvts_fl(curclass)
    [OUT.ts{i,j,curclass}(indkX), ...
        OUT.CVtargs{i,j,curclass}(:,indkX), ...
        OUT.CVdecs{i,j,curclass}(:,indkX), iModel] = nk_GetTestPerf(Ytest{curclass}{k(curclass)}, tCVL{curclass}, ...
                                                         [], iModel, Ymodel{curclass}{k(curclass)});
    if isnan(OUT.ts{i,j,curclass}(indkX))
        warning('Non-finite performance measures found in CV1 test data')
    end
    if VERBOSE, fprintf(', CV = %1.2f',OUT.ts{i,j,curclass}(indkX)); end
end

OUT.mdl{i,j,curclass}{indkX} = iModel;
% _________________________________________________________________________
function OUT = assign2out( OUT, RetrainImmediate, fFull, cvts_fl, ...
                            i , j, curclass)
    
%% Convert to single to save disk space
OUT.Trtargs{i,j,curclass}    = single(OUT.Trtargs{i,j,curclass}); 
OUT.Trdecs{i,j, curclass}    = single(OUT.Trdecs{i,j,curclass});
OUT.tr{i,j, curclass}        = single(OUT.tr{i,j,curclass}); 

% The decision scores for the CV1 data partition are only transferred to
% the OUT structure if models are NOT retrained OR if the CV1 training/test
% split is bypassed because optimization is not needed.
if ~fFull || RetrainImmediate || cvts_fl(curclass)
    OUT.CVtargs{i,j, curclass}   = single(OUT.CVtargs{i,j,curclass}); 
    OUT.CVdecs{i,j, curclass}    = single(OUT.CVdecs{i,j,curclass});
    OUT.ts{i,j, curclass}        = single(OUT.ts{i,j,curclass}); 
end
