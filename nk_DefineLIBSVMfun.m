function [LIBSVMTRAIN, LIBSVMPREDICT] = nk_DefineLIBSVMfun(SVM)

switch SVM.LIBSVM.LIBSVMver
    case 2
        LIBSVMTRAIN = 'svmtrain289';
        LIBSVMPREDICT = 'svmpredict289';
    case 0
        LIBSVMTRAIN = 'svmtrain312';
        LIBSVMPREDICT = 'svmpredict312';
    case 1
        LIBSVMTRAIN = 'svmtrain291';
        LIBSVMPREDICT = 'svmpredict291';
    case 3
        LIBSVMTRAIN = 'svmtrain289PLUS';
        LIBSVMPREDICT = 'svmpredict289PLUS';
end