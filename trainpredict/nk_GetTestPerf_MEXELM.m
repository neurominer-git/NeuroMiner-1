function [rs, ds] = nk_GetTestPerf_MEXELM(~, tXtest, ~, model, ~, ~)

scores = mexElmPredict( model.inW, model.bias, model.outW, tXtest');
ds = scores(1,:)'; rs = sign(ds); 