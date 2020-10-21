function [rs, ds] = nk_GetTestPerf_LIBLIN(~, tXtest, Ytest, md, ~, ~)
global SVM

[err_test, ~, predict_test] = predict_liblin22(Ytest, sparse(tXtest), md, sprintf(' -b %g -q',SVM.LIBLIN.b )); 
rs = err_test; ds = predict_test(:,1);

end