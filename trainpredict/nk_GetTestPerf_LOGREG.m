function [rs, ds] = nk_GetTestPerf_LOGREG(~, tXtest, ~, md, ~, ~)

[rs, ds] = LRpredict(md.theta, tXtest);
ds = ds(:,1);
