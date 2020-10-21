function algostr = getAlgoStr(analysis)

switch analysis.params.TrainParam.FUSION.flag
    case 3
        algostr = 'DESCFUSE';
    otherwise
        algostr = analysis.params.TrainParam.SVM.prog;
end
