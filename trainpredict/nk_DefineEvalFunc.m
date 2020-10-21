function funcstr = nk_DefineEvalFunc(GridParam)
global SVM

if ~exist('GridParam','var') || isempty(GridParam) && ~isemtpy(SVM.GridParam)
    GridParam = SVM.GridParam;
end
switch GridParam
    
    case 1  % Accuracy
        funcstr = 'ACCURACY';

    case 2  % True Positive Rate
        funcstr = 'TPR';

    case 3  % Sensitivity: TN / (TN + FP)
        funcstr = 'SENSITIVITY';

    case 4  % False Positive Rate = 1 - Specificity
        funcstr = 'FPR';

    case 5  % Positive Predictive Value: TP / (TP + FP)
        funcstr = 'PPV';

    case 6  % MCC
        funcstr = 'MCC';
      
    case 7 % AUC
        funcstr = 'AUC';
        
    case 8 % All parameters
        funcstr = 'ALLPARAM';
        
    case 9 % MSE
        funcstr = 'MSE';

    case 10 % Squared correlation coefficient
        funcstr = 'SCC';

    case 11 % Normalized root of MSE (NRMSD)
        funcstr = 'NRMSD';

    case 12 % Root of mean squared deviation (RMSD)
        funcstr = 'RMSD';
        
    case 13 % Gmean
        funcstr = 'GMEAN';
        
    case 14 % Balanced accuracy
        funcstr = 'BAC';
                
    case 15 % F-Score
        funcstr = 'FSCORE';
        
    case 16 % Correlation coefficient
        funcstr = 'CC';
        
    case 17 % Enhanced balanced accuracy (sens * spec) / 100
        funcstr = 'BAC2';
        
    case 18 % Mean average error
        funcstr = 'MAERR';
        
    case 19 % Prognostic Summary Index
        funcstr = 'PSI';
        
    case 20
        funcstr = 'NNP';
end

