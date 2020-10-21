function [rs, ds] = nk_GetTestPerf_MKLRVM(X, tXtest, Ytest, md, Features, k)
global MKLRVM

% Get Training Data
tX      = nk_ExtractFeatures(X, Features,[],k);

if iscell(tXtest), 
    funceval = [MKLRVM.funcname_predict '_MKL']; 
else
    funceval = MKLRVM.funcname_predict; 
end

t = nk_MakeRVMTargetMatrix(Ytest)';

err_prob = feval(funceval, md, tXtest, t, tX);

ds = err_prob(:,1); [dum,rs] = max(err_prob,[],2); rs(rs==2) =-1;    
        
end