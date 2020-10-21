function [mDTs, mTTs, Classes, mcolstart, mcolend] = nk_MultiAssemblePredictions( tsD, tsT, mDTs, mTTs, Classes, ul, curclass, mcolend )

ClassVec = ones(1,size(tsD,2))*curclass;

% Compute column pointers for multi-group CV2 array construction
mcolstart = mcolend + 1; mcolend = mcolstart + (ul-1);

% Enter decision values and predictions
% for multi-group classification into CV2 arrays
mDTs(:,mcolstart:mcolend) = tsD; 
mTTs(:,mcolstart:mcolend) = tsT;
Classes(1,mcolstart:mcolend) = ClassVec;
