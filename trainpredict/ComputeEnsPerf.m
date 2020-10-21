function perf = ComputeEnsPerf(perf, O)

nclass = size(O.TrHDperf,3);  
            
% Ensemble performance
perf.Ens_TrDecPerf = O.TrHDperf;
perf.Ens_CVDecPerf = O.CVHDperf;
perf.Ens_TrTargPerf = O.TrHTperf;
perf.Ens_CVTargPerf = O.CVHTperf;
perf.Ens_TrTargDiv = O.TrDiv;
perf.Ens_CVTargDiv = O.CVDiv;


    for curclass=1:nclass

        % Descriptive Stats
        TrDecPerf = reshape(O.TrHDperf(:,:,curclass),numel(O.TrHDperf(:,:,curclass)),1);
        TrTargPerf = reshape(O.TrHTperf(:,:,curclass),numel(O.TrHTperf(:,:,curclass)),1);
        CVDecPerf = reshape(O.CVHDperf(:,:,curclass),numel(O.CVHDperf(:,:,curclass)),1);
        CVTargPerf = reshape(O.CVHTperf(:,:,curclass),numel(O.CVHTperf(:,:,curclass)),1);
        TrTargDiv = reshape(O.TrDiv(:,:,curclass),numel(O.TrDiv(:,:,curclass)),1);
        CVTargDiv = reshape(O.CVDiv(:,:,curclass),numel(O.CVDiv(:,:,curclass)),1);

        % of CV1 training data
        perf.Ens_MeanTrDecPerf(curclass) = mean(TrDecPerf);
        perf.Ens_SDTrDecPerf(curclass) = std(TrDecPerf);
        perf.Ens_MeanTrDecPerf(curclass) = mean(TrTargPerf);
        perf.Ens_SDTrDecPerf(curclass) = std(TrTargPerf);
        perf.Ens_MeanTrDiv(curclass) = mean(TrTargDiv);
        perf.Ens_SDTrDiv(curclass) = std(TrTargDiv);

        % of CV1 test data
        perf.Ens_MeanCVDecPerf(curclass) = mean(CVDecPerf);
        perf.Ens_SDCVDecPerf(curclass) = std(CVDecPerf);
        perf.Ens_MeanCVTargPerf(curclass) = mean(CVTargPerf);
        perf.Ens_SDCVTargPerf(curclass) = std(CVTargPerf);
        perf.Ens_MeanCVDiv(curclass) = mean(CVTargDiv);
        perf.Ens_SDCVDiv(curclass) = std(CVTargDiv);
    end

end
