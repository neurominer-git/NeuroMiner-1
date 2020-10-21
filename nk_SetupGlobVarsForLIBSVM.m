function [SVM, CMDSTR, LIBSVMTRAIN, LIBSVMPREDICT] = nk_SetupGlobVarsForLIBSVM(param)
global MODEFL

SVM.prog = 'LIBSVM';
SVM.LIBSVM = param.LIBSVM;
SVM.kernel = param.kernel;
SVM.GridParam = param.evalfunc;
CMDSTR = nk_DefineCmdStr(SVM, MODEFL);

switch SVM.LIBSVM.LIBSVMver
    case 0
        LIBSVMTRAIN = 'svmtrain312';
        LIBSVMPREDICT = 'svmpredict312';
    case 2
        LIBSVMTRAIN = 'svmtrain289';
        LIBSVMPREDICT = 'svmpredict289';
    case 1
        LIBSVMTRAIN = 'svmtrain291';
        LIBSVMPREDICT = 'svmpredict291';
    case 3
        LIBSVMTRAIN = 'svmtrain289PLUS';
        LIBSVMPREDICT = 'svmpredict289PLUS';
end