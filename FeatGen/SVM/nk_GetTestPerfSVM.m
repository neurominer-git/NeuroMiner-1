function [rs, ds] = nk_GetTestPerfSVM(tXtest, Ytest, md, SVM, LIBSVMPREDICT)

[err_test, dum, predict_test] = feval( LIBSVMPREDICT, Ytest, tXtest, md, sprintf(' -b %g',SVM.LIBSVM.Optimization.b)); 
rs = err_test; ds = predict_test(:,1);

end