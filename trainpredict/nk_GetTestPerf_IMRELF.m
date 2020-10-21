function [rs, ds] = nk_GetTestPerf_IMRELF(~, tXtest, ~, md, ~, ~)

[err_test, predict_test] = ...
    IMRelief_Sigmoid_FastImple_Predict2(md.X, md.targets, tXtest', md.Weight, ...
                                        md.distance, md.sigma);
rs = err_test; ds = predict_test(:,1);

end