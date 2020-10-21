function [entrop, perf,D, T, L] = nk_ComputeEntropy(decs, targs, TsInd, TsIdx, TsLabels)
global EVALFUNC

[D, T, L] = nk_Vals2Ind(decs, targs, TsInd, TsIdx, TsLabels);

lbind = any(T,2);
hx = sign(sum(T(lbind,:),2));
lb = sign(sum(L(lbind,:),2));
perf = feval(EVALFUNC,hx,lb);
entrop = nk_Ambiguity(T(lbind,:)); 

return