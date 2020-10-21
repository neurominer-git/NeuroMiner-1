function [ Tcorr, Dcorr ] = nk_ApplyDetrend(XTest, YTrain, XTrain, Mkl, Fkl, detrend)

global MODEFL

nMkl = numel(Mkl); nY = numel(YTrain);
Tcorr = zeros(nY, nMkl); Dcorr = zeros(nY, nMkl); 
for z = 1 : nMkl     
    [~, Tcorr(:,z), Dcorr(:,z)] = nk_GetTestPerf(XTest, YTrain, Fkl(:,z), Mkl{z}, XTrain, 1);
    switch MODEFL
        case 'regression'
            Tcorr(:,z) = nk_DetrendPredictions2(detrend.beta, detrend.p, Tcorr(:,z)); Dcorr(:,z) = Tcorr(:,z);
        case 'classification'
            Dcorr(:,z) = Dcorr(:,z) - detrend.thresh; Tcorr(:,z) = sign(Dcorr(:,z));
    end
            
end
