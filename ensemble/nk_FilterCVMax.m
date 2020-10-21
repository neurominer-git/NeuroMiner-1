function [Fmask, opt_heD, opt_heT, opt_D, opt_T] = nk_FilterCVMax(D, T, L, F, kSub, kSp, kInd, MinNum)

% Select classifiers acording to their decision values / target predictions
% using the margin-based g-flip algorithm

% First, 
[heD, heT] = nk_EnsPerf(D, T, L);


return