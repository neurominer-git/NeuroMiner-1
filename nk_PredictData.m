% =========================================================================
% function Predict = nk_PredictData(F, W, TR, TRInd, dTRLabel, ...
%                                       CVD, CVDInd, dCVDLabel, ...
%                                       TS, dTSInd, dTSLabel, ...
%                                       mTSInd, mTSLabel, ...
%                                       MD, detrend)
% =========================================================================
% This function applies a CV1 ensemble of models to the CV2 validation
% data.
%
% INPUTS:
% -------
% 
% OUTPUTS:
% --------
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Nikolaos Koutsouleris, 07/2017

function Predict = nk_PredictData(F, W, TR, TRInd, dTRLabel, ...
                                        CVD, CVDInd, dCVDLabel, ...
                                        TS, dTSInd, dTSLabel, ...
                                        mTSInd, mTSLabel, ...
                                        MD, ngroups, detrend)
                                    
global SVM CV MULTI MODEFL VERBOSE RFE EVALFUNC RAND MULTILABEL

% Initialize Variables
[iy,jy]     = size(TR);
dTRL    = dTRLabel{1,1}{1}(:,MULTILABEL.curdim);
lgroups = unique(dTRL(~isnan(dTRL)));

if strcmp(MODEFL,'classification')
    if RAND.Decompose ~= 9
        nclass = length(CV.class{1,1});
    else
        nclass = 1;
    end
elseif strcmp(MODEFL,'regression')
    nclass = 1;
end
nF          = size(F,3);
nW          = size(W,3);
ENSMETHOD   = 1;

% CV1 partitions
Predict.binCV1Predictions        = cell(iy,jy);
Predict.binCV1Performance        = zeros(iy,jy,nclass);
Predict.binCV1Performance_Mean   = zeros(1,nclass);
Predict.binCV1Performance_SD     = zeros(1,nclass);
if MULTI.flag, 
    Predict.MultiCV1Predictions  = cell(iy,jy);
    Predict.MultiCV1Probabilities  = cell(iy,jy);
    Predict.MultiCV1Performance  = zeros(iy,jy); 
end

switch SVM.prog
    case 'SEQOPT'
        Predict.binCV1CasePropagations = cell(iy,jy,nclass);
        Predict.binCV1PerformanceIncreases = cell(iy,jy,nclass);
        Predict.binCV1DecValTraj = cell(iy,jy,nclass);
    case 'WBLCOX'
        Predict.binCV1times = cell(iy,jy,nclass);
        Predict.binCV1probthresh = zeros(iy,jy,nclass);
end

% Compute size of multi-group and binary arrays to avoid dynamic memory allocation
% This significantly improves code execution performance
Ydims = nk_GetCV2EnsembleDims(F); if nF ~= nclass, sYdims = Ydims*nclass; else, sYdims = sum(Ydims); end
if iscell(TS{1,1})
    Xdims = size(TS{1,1}{1},1);
else
    Xdims = size(TS{1,1},1);
end
mDTs = zeros(Xdims, sYdims); mTTs = mDTs; Classes = zeros(1,sYdims);
mcolend = 0;

% Check if detrending should be applied after label prediction
if exist('detrend','var') && ~isempty(detrend)
    detrendfl = true;
else
    detrendfl = false;
end

for k=1:iy % Loop through CV1 permutations

    for l=1:jy % Loop through CV1 folds
        
        for curclass = 1:nclass % Loop through dichotomizers
            
            % Extract Test Data
            if iscell(TS{k,l}), 
                if ~isempty(mTSInd),
                    if iscell(mTSInd{k,l})
                        XTest = TS{k,l}{curclass}(mTSInd{k,l}{curclass},:); 
                    else
                        XTest = TS{k,l}{curclass}(mTSInd{k,l},:); 
                    end
                else
                    XTest = TS{k,l}{curclass}; 
                end
            else
                if ~isempty(mTSInd)
                    if iscell(mTSInd{k,l})
                        XTest = TS{k,l}(mTSInd{k,l}{curclass},:);
                    else
                        XTest = TS{k,l}(mTSInd{k,l},:);
                    end
                else
                    XTest = TS{k,l};
                end
            end
            
            % Extract Training Data
            if iscell(TR{k,l}), 
                if ~isempty(TRInd), 
                    if iscell(TRInd{k,l})
                        XTrain = TR{k,l}{curclass}(TRInd{k,l}{curclass},:); 
                    else
                        XTrain = TR{k,l}{curclass}(TRInd{k,l},:); 
                    end
                else
                    XTrain = TR{k,l}{curclass}; 
                end
            else
                if ~isempty(TRInd)
                    if iscell(TRInd{k,l})
                        XTrain = TR{k,l}(TRInd{k,l}{curclass},:);
                    else
                        XTrain = TR{k,l}(TRInd{k,l},:);
                    end
                else
                    XTrain = TR{k,l};
                end
            end
            YTrain = dTRLabel{k,l}{curclass}(:,MULTILABEL.curdim);
            
            %Extract CV1 Test Data
            if iscell(CVD{k,l}), 
                if ~isempty(CVDInd), 
                    if iscell(CVDInd{k,l})
                        XCV = CVD{k,l}{curclass}(CVDInd{k,l}{curclass},:); 
                    else
                        XCV = CVD{k,l}{curclass}(CVDInd{k,l},:); 
                    end
                else
                    XCV = CVD{k,l}{curclass}; 
                end
            else
                if ~isempty(CVDInd)
                    if iscell(CVDInd{k,l})
                        XCV = CVD{k,l}(CVDInd{k,l}{curclass},:);
                    else
                        XCV = CVD{k,l}(CVDInd{k,l},:);
                    end
                else
                    XCV = CVD{k,l};
                end
            end
            YCV = dCVDLabel{k,l}{curclass}(:,MULTILABEL.curdim);
            Ydum = zeros(size(XTest,1),1);
            if RFE.ClassRetrain
               XTrain = [XTrain; XCV]; YTrain = [YTrain; YCV]; 
            end
            
            if nF == 1
                ul = size(F{k,l},2); Fkl = F{k,l};
            else
                ul = size(F{k,l,curclass},2); Fkl = F{k,l,curclass}; 
            end
            
            %%%%%%%%%%%%%%%% GET FEATURE SUBSPACE MASK %%%%%%%%%%%%%%%%
            if ~islogical(Fkl), Fkl = Fkl ~= 0; end
            
            % Optionally, retrain model using XTrain if no predefined model
            % exists
            Mkl = MD{k,l,curclass};
            
            %%%%%%%%%%%%%%%%%%% PREDICT TEST DATA %%%%%%%%%%%%%%%%%%%%%
            if VERBOSE, 
                switch MODEFL
                    case 'classification'
                        if RAND.Decompose ~= 9
                            fprintf('\nPredicting data in CV1 [%g,%g], classifier #%g (%s), %g models.', ...
                                k, l, curclass, CV.class{1,1}{curclass}.groupdesc, ul); 
                        else
                            fprintf('\nPredicting data in CV1 [%g,%g], multi-group classifier, %g models.', ...
                                k, l, ul); 
                        end
                    case 'regression'
                        fprintf('\nPredicting data in CV1 [%g,%g], regressor, %g models.', ...
                        k, l, ul); 
                end
            end
            
            if detrendfl
               [ tsT, tsD ] = nk_ApplyDetrend(XTest, Ydum, XTrain, Mkl, Fkl, detrend);
            else
                [~, tsT, tsD, Mkl] = nk_GetTestPerf(XTest, Ydum, Fkl, Mkl, XTrain, 1);
            end
               
            Predict.binCV1Predictions{k,l,curclass} = tsD;
            switch SVM.prog
                case 'SEQOPT'
                    for u=1:size(tsD,2)
                        Predict.binCV1CasePropagations{k,l,curclass} = [Predict.binCV1CasePropagations{k,l,curclass} Mkl{u}.Nremain_test];
                        Predict.binCV1PerformanceIncreases{k,l,curclass} = [Predict.binCV1PerformanceIncreases{k,l,curclass}; Mkl{u}.SeqPerfGain_test];
                    end
                    Predict.binCV1DecValTraj{k,l,curclass} = Mkl{u}.optDh;
                case 'WBLCOX'
                    mdthresh = zeros(size(tsD,2),1);
                    for u=1:size(tsD,2)
                        Predict.binCV1times{k,l,curclass} = [ Predict.binCV1times{k,l,curclass} Mkl{u}.predicted_time]; 
                        mdthresh(u) = Mkl{u}.cutoff.prob;
                    end
                    Predict.binCV1probthresh(k,l,curclass) = nm_nanmedian(mdthresh(u));
            end
            n_subj = size(tsD,1);
            
            % Weight decision values
            if ~isempty(W) && ~isempty(W{k,l,curclass}), 
                if nW == 1 , Wkl = W{k,l}';  else Wkl = W{k,l,curclass}'; end
                wx = repmat(Wkl,n_subj,1); tsD = bsxfun(@times,tsD,wx); tsT = bsxfun(@times,tsT,wx); 
            end
            
            % Compute binary performance on CV2 test data for current CV1 partition
            if ~isempty(tsD) 
                dtsD = tsD(dTSInd{curclass},:); 
            else
                tsD = 0; tsT = 0; dtsD = [];
            end
            perf = zeros(ul,1);
            
            for u=1:ul
                if isempty(dtsD)
                    perf(u) = 0;
                else
                    offs =0; if strcmp(SVM.prog,'WBLCOX'), offs =  Mkl{u}.cutoff.prob; end
                    perf(u) = feval(EVALFUNC, dTSLabel{curclass}(:,MULTILABEL.curdim), dtsD(:,u)-offs);
                end
            end
            
            Predict.binCV1Performance(k,l,curclass) = mean(perf);

            % Multi-group CV2 array construction
            [mDTs, mTTs, Classes, mcolstart, mcolend] = ...
                nk_MultiAssemblePredictions( tsD, tsT, mDTs, mTTs, Classes, ul, curclass, mcolend );
            
            if curclass == 1, mcolX = mcolstart; end
            
        end
        
        % Compute multi-group performance on CV2 test data for current CV1 partition
        if MULTI.flag
            [Predict.MultiCV1Performance(k,l), Predict.MultiCV1Predictions{k,l}, Predict.MultiCV1Probabilities{k,l}] = ...
                nk_MultiEnsPerf(mDTs(:,mcolX:mcolend), ...
                mTTs(:,mcolX:mcolend), ...
                mTSLabel(:,MULTILABEL.curdim), ...
                Classes(:,mcolX:mcolend), ngroups);
        end
    end
end

% ScaleFlag = false;

% Compute CV2 ensemble binary performance for current CV2 partition
for curclass = 1 : nclass
    
    % Extract current dichotomization decision from (multi-group) array
    indCurClass = Classes == curclass;
    if isempty(dTSInd{curclass}), 
        Predict.binCV2Performance_Targets(curclass) = NaN;
        Predict.binCV2Performance_DecValues(curclass) = NaN;   
        Predict.binCV2Diversity_Targets(curclass) = NaN;
        Predict.binCV2Diversity_DecValues(curclass) = NaN;
        Predict.binCV2Predictions_Targets{curclass} = NaN;
        Predict.binCV2Predictions_DecValues{curclass} = NaN;
        continue; 
    end 
    dDTs = mDTs(dTSInd{curclass},indCurClass);
    dTTs = mTTs(dTSInd{curclass},indCurClass);
    perf = Predict.binCV1Performance(:,:,curclass);
    Predict.binCV1Performance_Mean(1,curclass) = nm_nanmean(perf(:));
    Predict.binCV1Performance_SD(1,curclass) = nm_nanstd(perf(:));

    switch MODEFL

        case 'classification'
            
            switch RAND.Decompose
                
                case 9 % Multi-group classifier
                    
                    % Generate majority voting arry
                    Tvotes = zeros(size(dTTs,1),ngroups);
                    
                    for i = 1:ngroups
                       
                        Tvotes(:,i) = nm_nansum(dTTs == lgroups(i),2);
                        
                    end
                    [~,ind] = max(Tvotes,[],2);
                    hrx = lgroups(ind);
                    hdx = hrx;
                    
                otherwise
            
                    switch ENSMETHOD

                        case 1 % This is the simple majority vote
                            % Hard label decision:
                            hrx = sign(nm_nansum(dTTs,2));
                            hdx = nm_nanmedian(dDTs,2);

                        case 2 % Product majority vote
                            hrx = sign(prod(dTTs,2));
                            hdx = prod(dDTs,2);

                        case 3 % Error Correcting Output Codes
                            coding=1; decoding=1;
                            classes = ones(1,size(dTTs,2));
                            %hrs(hrs==-1) = ;
                            hrx = nk_ErrorCorrOutCodes(dTTs, classes, coding, decoding);
                            hrx(hrx==2) = -1; hdx = hrx;
                    end
            end
            % Check for zeros
            if sum(hrx==0) > 0 % Throw coin
                hrx = nk_ThrowCoin(hrx);
            end 
            
            Predict.binCV2Performance_Targets(curclass) = feval(EVALFUNC, dTSLabel{curclass}(:,MULTILABEL.curdim), hrx);
            mdcutoff=0; if strcmp(SVM.prog,'WBLCOX'), probthresh = Predict.binCV1probthresh(:,:,curclass); mdcutoff = nm_nanmedian(probthresh(:)); end
            Predict.binCV2Performance_DecValues(curclass) = feval(EVALFUNC, dTSLabel{curclass}(:,MULTILABEL.curdim), hdx-mdcutoff);
            % TODO: Adapt the next coding statement to feval
            Predict.binCV2Diversity_Targets(curclass) = nk_Entropy(dTTs);
            Predict.binCV2Diversity_DecValues(curclass) = nk_RegAmbig(dDTs);

        case 'regression'
            hrx     = nm_nanmedian(dTTs,2);
            hdx     = hrx;
            
            Predict.binCV2Performance_Targets(curclass) = feval(EVALFUNC, dTSLabel{curclass}(:,MULTILABEL.curdim), hrx);
            Predict.binCV2Performance_DecValues(curclass) = Predict.binCV2Performance_Targets(curclass);
            Predict.binCV2Diversity_Targets(curclass) = 0;
            Predict.binCV2Diversity_DecValues(curclass) = 0;
            
    end
    
    Predict.binCV2Predictions_Targets{curclass} = hrx;
    Predict.binCV2Predictions_DecValues{curclass} = hdx;
    
end

% Compute CV2 ensemble multi-group performance and class membership
if MULTI.flag  
    [Predict.MultiCV2Performance, Predict.MultiCV2Predictions, Predict.MultiCV2Probabilities] = ...
                nk_MultiEnsPerf(mDTs, mTTs, mTSLabel(:,MULTILABEL.curdim), Classes, ngroups );
    Predict.MultiCV2Diversity_Targets = nm_nanmean(Predict.binCV2Diversity_Targets);
    Predict.MultiCV2Diversity_DecValues = nm_nanmean(Predict.binCV2Diversity_DecValues);
end

return