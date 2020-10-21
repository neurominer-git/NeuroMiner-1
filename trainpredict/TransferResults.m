%==========================================================================
%FORMAT o = TransferResults( I, O, Param )                                    
%==========================================================================
%Subfunction of OptimCore
%
%Generates output array from OptimCore processing results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%(c) Nikolaos Koutsouleris, 03 / 2020
function o = TransferResults(I, O, Param)

global MULTI W2AVAIL 
                                
o.SubSpaces                 = I.F;
o.Models                    = O.mdl;
if Param.SubSpaceFlag
    o.SubSpacesMask             = O.F;
    o.Weights                   = O.Weights;
    o.TrPredictions             = O.TrHT;
    o.TrDecisionValues          = O.TrHD;
    o.CVPredictions             = O.CVHT;
    o.CVDecisionValues          = O.CVHD;
    o.TrHTperf                  = O.TrHTperf;
    o.TrHDperf                  = O.TrHDperf;
    o.CVHTperf                  = O.CVHTperf;
    o.CVHDperf                  = O.CVHDperf;
else
    o.SubSpaceMask              = O.F;
    o.TrPredictions             = O.Trtargs;
    o.TrDecisionValues          = O.Trdecs;
    o.CVPredictions             = O.CVtargs;
    o.CVDecisionValues          = O.CVdecs;
    o.TrHTperf                  = cell2mat(O.tr);
    o.TrHDperf                  = cell2mat(O.tr);
    o.CVHTperf                  = cell2mat(O.ts);
    o.CVHDperf                  = cell2mat(O.ts);
end

if W2AVAIL
    o.w2 = O.w2;
    o.Md = O.Md;
    o.Mm = O.Mm;
end

for curclass = 1 : I.nclass
    [ix, jx] = size(o.TrHTperf(:,:,curclass));
    o.MeanTrHTperf(curclass) = nm_nanmean( reshape(o.TrHTperf(:,:,curclass),ix*jx,1) );
    o.SDTrHTperf(curclass)   = nm_nanstd( reshape(o.TrHTperf(:,:,curclass),ix*jx,1) );
    o.MeanTrHDperf(curclass) = nm_nanmean( reshape(o.TrHDperf(:,:,curclass),ix*jx,1) );
    o.SDTrHDperf(curclass)   = nm_nanstd( reshape(o.TrHDperf(:,:,curclass),ix*jx,1) );
    
    o.MeanCVHTperf(curclass) = nm_nanmean( reshape(o.CVHTperf(:,:,curclass),ix*jx,1) );
    o.SDCVHTperf(curclass)   = nm_nanstd( reshape(o.CVHTperf(:,:,curclass),ix*jx,1) );
    o.MeanCVHDperf(curclass) = nm_nanmean( reshape(o.CVHDperf(:,:,curclass),ix*jx,1) );
    o.SDCVHDperf(curclass)   = nm_nanstd( reshape(o.CVHDperf(:,:,curclass),ix*jx,1) );
    
    if Param.SubSpaceFlag
        o.Ens_MeanTrDiv(curclass) = nm_nanmean(reshape(O.TrDiv(:,:,curclass),numel(O.TrDiv(:,:,curclass)),1));
        o.Ens_MeanCVDiv(curclass) = nm_nanmean(reshape(O.CVDiv(:,:,curclass),numel(O.CVDiv(:,:,curclass)),1));
    
        o.Ens_SDTrDiv(curclass) = nm_nanstd(reshape(O.TrDiv(:,:,curclass),numel(O.TrDiv(:,:,curclass)),1));
        o.Ens_SDCVDiv(curclass) = nm_nanstd(reshape(O.CVDiv(:,:,curclass),numel(O.CVDiv(:,:,curclass)),1));
    end
end

if MULTI.flag
    if Param.SubSpaceFlag
        o.MultiTrPredictions    = O.mTrPred;
        o.MultiCVPredictions    = O.mCVPred;          
        o.MultiTrPerf           = O.mTrPerf;
        o.MultiCVPerf           = O.mCVPerf;
        o.Ens_MeanMultiTrDiv    = mean(O.mTrDiv(:));
        o.Ens_SDMultiTrDiv      = std(O.mTrDiv(:));
        o.Ens_MeanMultiCVDiv    = mean(O.mCVDiv(:));
        o.Ens_SDMultiCVDiv      = std(O.mCVDiv(:));
    else
        o.MultiTrPredictions    = O.mTrPred;
        o.MultiCVPredictions    = O.mCVPred;          
        o.MultiTrPerf           = cell2mat(O.mtr);
        o.MultiCVPerf           = cell2mat(O.mts);
    end
    o.MeanMultiTrPerf       = mean(o.MultiTrPerf(:));
    o.SDMultiTrPerf         = std(o.MultiTrPerf(:));
    o.MeanMultiCVPerf       = mean(o.MultiCVPerf(:));
    o.SDMultiCVPerf         = std(o.MultiCVPerf(:));
end

if isfield(O,'detrend') && O.detrend.flag
   o.detrend = O.detrend;
end

if isfield(O,'critgain')
    td = cell2mat(reshape(O.examfreq(:,:,1),ix*jx,1));
    o.MeanCritGain = zeros(I.nclass, size(td,2)-1);
    o.SDCritGain = zeros(I.nclass, size(td,2)-1);
    o.MeanExamFreq = zeros(I.nclass, size(td,2));
    o.SDExamFreq = zeros(I.nclass, size(td,2));
    o.MeanAbsThreshU = zeros(I.nclass, size(td,2)-1);
    o.SDAbsThreshU = zeros(I.nclass, size(td,2)-1);
    o.MeanAbsThreshL = zeros(I.nclass, size(td,2)-1);
    o.SDAbsThreshL = zeros(I.nclass, size(td,2)-1);
    o.MeanPercThreshU =  zeros(I.nclass, size(td,2)-1);
    o.SDPercThreshU =  zeros(I.nclass, size(td,2)-1);
    o.MeanPercThreshL =  zeros(I.nclass, size(td,2)-1);
    o.SDPercThreshL =  zeros(I.nclass, size(td,2)-1);
    for curclass = 1 : I.nclass
       o.CritGain = O.critgain;
       o.MeanCritGain(curclass,:) = nm_nanmean(cell2mat(reshape(O.critgain(:,:,curclass),ix*jx,1)));
       o.SDCritGain(curclass,:) = nm_nanstd(cell2mat(reshape(O.critgain(:,:,curclass),ix*jx,1)));
       o.MeanExamFreq(curclass,:) = nm_nanmean(cell2mat(reshape(O.examfreq(:,:,curclass),ix*jx,1)));
       o.SDExamFreq(curclass,:) = nm_nanstd(cell2mat(reshape(O.examfreq(:,:,curclass),ix*jx,1)));
       o.MeanAbsThreshU(curclass,:) = nm_nanmean(cell2mat(reshape(O.absthreshU(:,:,curclass),ix*jx,1)));
       o.SDAbsThreshU(curclass,:) = nm_nanstd(cell2mat(reshape(O.absthreshU(:,:,curclass),ix*jx,1)));
       o.MeanAbsThreshL(curclass,:) = nm_nanmean(cell2mat(reshape(O.absthreshL(:,:,curclass),ix*jx,1)));
       o.SDAbsThreshL(curclass,:) = nm_nanstd(cell2mat(reshape(O.absthreshL(:,:,curclass),ix*jx,1)));
       o.MeanPercThreshU(curclass,:) = nm_nanmean(cell2mat(reshape(O.percthreshU(:,:,curclass),ix*jx,1)));
       o.SDPercThreshU(curclass,:) = nm_nanstd(cell2mat(reshape(O.percthreshU(:,:,curclass),ix*jx,1)));
       o.MeanPercThreshL(curclass,:) = nm_nanmean(cell2mat(reshape(O.percthreshL(:,:,curclass),ix*jx,1)));
       o.SDPercThreshL(curclass,:) = nm_nanstd(cell2mat(reshape(O.percthreshL(:,:,curclass),ix*jx,1)));
    end
elseif isfield(O,'threshprob')
    td = cell2mat(reshape(O.times(:,:,1),ix*jx,1));
    o.PredictedTimes = O.times;
    o.MeanThreshProb = zeros(I.nclass,1);
    o.SDThreshProb = zeros(I.nclass,1);
    o.MeanThreshPerc = zeros(I.nclass,1);
    o.SDThreshPerc = zeros(I.nclass,1);
    o.MeanPredictedTimes = zeros(I.nclass, size(td,2));
    o.SDPredictedTimes = zeros(I.nclass, size(td,2));
    for curclass = 1 : I.nclass
        o.MeanThreshProb(curclass) = nm_nanmean(reshape(O.threshprob(:,:,curclass),ix*jx,1));
        o.SDThreshProb(curclass) = nm_nanstd(reshape(O.threshprob(:,:,curclass),ix*jx,1));
        o.MeanThreshPerc(curclass) = nm_nanmean(reshape(O.threshperc(:,:,curclass),ix*jx,1));
        o.SDThreshPerc(curclass) = nm_nanstd(reshape(O.threshperc(:,:,curclass),ix*jx,1));
        o.MeanPredictedTimes(curclass,:) = nm_nanmean(cell2mat(reshape(O.times(:,:,curclass),ix*jx,1)));
        o.SDPredictedTimes(curclass,:) = nm_nanstd(cell2mat(reshape(O.times(:,:,curclass),ix*jx,1)));
    end
end

