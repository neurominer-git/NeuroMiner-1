function [trainfunc, predictfunc] = nk_DefineTrainPredictFunc(experimental)

global SVM

if exist('experimental','var') && ~isempty(experimental) && experimental
    TrainBase = 'nk_GetParam2_';
    PredictBase = 'nk_GetTestPerf_';
else
    TrainBase = 'nk_GetParam_';
    PredictBase = 'nk_GetTestPerf_';
end

trainfunc = [ TrainBase SVM.prog ];
predictfunc = [ PredictBase SVM.prog ] ;
      
end