function GD = nk_CompPlattProbOptimum(GD, Ps, mapY)
global CV SVM MULTI RFE

%%%% PREPARATIONS %%%%
[CV1perms, CV1folds] = size(mapY.Tr);
nclass = length(CV.class{1,1});
F = cell(CV1perms, CV1folds, nclass); W = F; 
SlackParam = cell(nclass,1); KernParam = cell(nclass,1);

%%%% Loop through dichotomizers %%%%
for curclass=1:nclass
    
    % THIS HAS TO BE REWORKED!
%     cPs = cell(nclass,1);
%     i = GD.BinaryGridSelection.
%     for curclass = 1:nclass
% 
%         cPs{curclass} = Ps{curclass}(i,:); 
%         
%         switch SVM.prog
%              case {'MKLRVM'}
%                  if nvar >1
%                     cPs{curclass} = repmat({cPs{curclass}},1,nvar);
%                  end
%             case {'LIBSVM','LIBLIN','SVMLIT'}
%                  % Convert parameters to char array
%                  cPs{curclass} = num2str(cPs{curclass}','%1.10f');
%                  % Concatenate parameter string
%                  cPs{curclass} = nk_ConcatLIBSVMParamStr(cPs{curclass});
%             case 'IMRELF'
% 
%         end
% 
%     end
    
    nc = size(GD.BinaryGridSelection{curclass}.bestfeats,3);
    if nc > 1
        F(:,:,curclass) = GD.BinaryGridSelection{curclass}.bestfeats(:,:,curclass);
        W(:,:,curclass) = GD.BinaryGridSelection{curclass}.bestweights(:,:,curclass);
    else
        F(:,:,curclass) = GD.BinaryGridSelection{curclass}.bestfeats{1};
        W(:,:,curclass) = GD.BinaryGridSelection{curclass}.bestweights{1};
    end
end

TR = mapY.Tr;   TRInd = mapY.TrInd; dTRLabel = mapY.TrL; 
CVD = mapY.CV;  CVDInd = mapY.CVInd; dCVDLabel = mapY.CVL;
TS = mapY.Ts;   dTSInd = mapY.TsInd; mTSInd = []; dTSLabel = mapY.TsL; 
if ~isempty(MULTI) && MULTI.flag, mTSLabel = mapY.mTsL; else mTSLabel = []; end
MD = []; 

%%%% Recompute CV2 predictions using probability flag = true %%%%
BinCV2results = nk_PredictData(F, W, TR, TRInd, dTRLabel, ...
                                        CVD, CVDInd, dCVDLabel, ...
                                        TS, dTSInd, dTSLabel, ...
                                        mTSInd, mTSLabel, ...
                                        MD, 0, SlackParam, KernParam, true);

            
for curclass=1:nclass
    GD.BinaryGridSelection{curclass}.bestpred = BinCV2results.binCV1Predictions(:,:,curclass);
    switch RFE.CV2Class.type 
        case 1
            GD.BinaryGridSelection{curclass}.besttestparam = BinCV2results.BinCV1Performance_Mean(1,curclass);
        case 2
            switch RFE.CV2Class.EnsembleStrategy.Metric                                     
                case 1
                    GD.BinaryGridSelection{curclass}.besttestparam = BinCV2results.binCV2Performance_Targets(curclass);
                case 2
                    GD.BinaryGridSelection{curclass}.besttestparam = BinCV2results.binCV2Performance_DecValues(curclass);
            end
    end
end
                                    
                                    
if ~isempty(MULTI) && MULTI.flag 
    
    if ~MULTI.BinBind
    
        SlackParam = 0; KernParam = 0;

        if isfield(GD.MultiGroupGridSelection,'bestc')
            SlackParam = GD.MultiGroupGridSelection.bestc;
        end

        if isfield(GD.MultiGroupGridSelection,'bestg')
            KernParam = GD.MultiGroupGridSelection.bestg;        
        end

        switch SVM.prog
            case {'CUDSVM','MikRVM','MKLRVM','GLMFIT','LIKNON', 'kNNMEX'}
            otherwise
                SlackParam = cellfun(@(x) num2str(x,'%1.10f'),SlackParam,'UniformOutput',false);
                KernParam = cellfun(@(x) num2str(x,'%1.10f'),KernParam,'UniformOutput',false);       
        end

        F = GD.MultiGroupGridSelection.bestfeats;
        W = GD.MultiGroupGridSelection.bestweights;

        MultiCV2results = nk_PredictData(F, W, TR, TRInd, dTRLabel, ...
                                            CVD, CVDInd, dCVDLabel, ...
                                            TS, dTSInd, dTSLabel, ...
                                            mTSInd, mTSLabel, ...
                                            MD, repmat(SlackParam,nclass,1), repmat(KernParam,2,1), true);

        GD.MultiGroupGridSelection.bestpred = MultiCV2results.MultiCV2Predictions;
        switch RFE.CV2Class.type 
            case 1
                GD.MultiGroupGridSelection.besttestparam = mean(MultiCV2results.MultiCV1Performance(:));
            case 2
                GD.MultiGroupGridSelection.besttestparam = MultiCV2results.MultiCV2Performance;
        end
    else
        GD.MultiGroupGridSelection.bestpred = BinCV2results.MultiCV2Predictions;
        switch RFE.CV2Class.type 
            case 1
                GD.MultiGroupGridSelection.besttestparam = mean(BinCV2results.MultiCV1Performance(:));
            case 2
                GD.MultiGroupGridSelection.besttestparam = BinCV2results.MultiCV2Performance;
        end
    end
    
end
    
return

