function [rs, ds] = nk_GetTestPerf_ROBSVM(~, tXtest, Ytest, md, ~, ~)
global SVM 

if isfield(SVM.ROBSVM,'Optimization'), b= SVM.ROBSVM.Optimization.b; else, b=0; end
[err_test, ~, predict_test] = svmpredict312( Ytest, tXtest, md, sprintf(' -b %g', b)); 
rs = err_test; ds = predict_test(:,1);

