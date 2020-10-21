function [L, X, IND] = nk_LabelImputer( L, X, IND, Params, IMPUTE )
global PREDICTFUNC

if IMPUTE.flag 
    IMPUTE.X = X;
    for i=1:size(L,2) % Loop through labels
        if strcmp(IMPUTE.method,'ml')
            % Do ML-based label propagation by means of the user-defined
            % algo
            indf = ~isnan(L(:,i));
            X_tr = X(indf,:); L_tr = L(indf,i); X_ts = X(~indf,:); L_ts = zeros(sum(~indf),1);
            [~, model] = nk_GetParam2( X_tr, L_tr, Params, 1); 
            L(~indf,i)= feval(PREDICTFUNC, X_tr, X_ts, L_ts, model);
        elseif ~strcmp(IMPUTE.method,'none')
            % Do NN-based label propagation
            L(:,i) = nk_PerfImputeLabelObj(L(:,i), IMPUTE); 
        end
    end
else
    indf = sum(isnan(L),2)>0;
    L(indf,:)=[];
    X(indf,:)=[];
    IND(indf) = [];
end
