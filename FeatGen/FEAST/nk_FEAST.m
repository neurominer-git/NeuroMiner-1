function [ind, R] = nk_FEAST(Y, label, k, params)

global MODEFL VERBOSE

[~,n] = size(Y);

switch params.cmd

    case 'discretize'
        
        if VERBOSE; fprintf(' Discretize [ %g : %g : %g ]... ', ...
                params.DISCRET.binstart,params.DISCRET.binsteps, params.DISCRET.binstop); 
        end
        
        Y = discretize(Y, ...
            params.DISCRET.binstart, ...
            params.DISCRET.binsteps, ...
            params.DISCRET.binstop);
        
        if strcmp(MODEFL,'regression'), label = discretize(label, ...
            params.DISCRET.binstart, ...
            params.DISCRET.binsteps, ...
            params.DISCRET.binstop); 
        end
        
    case 'symbolize'
        
        if VERBOSE; fprintf(' Symbolize [ Min=%g, Max=%g, L=%g, #STD=%g ]... ', ...
                params.SYMBOL.symMinBin, params.SYMBOL.symMaxBin, params.SYMBOL.symSeqLength, params.SYMBOL.symStdNum); 
        end
        
        Y = symbolize(Y, params.SYMBOL.symMinBin, params.SYMBOL.symMaxBin, params.SYMBOL.symSeqLength, params.SYMBOL.symStdNum);
        
        if strcmp(MODEFL,'regression'), label = symbolize(label, ...
                params.SYMBOL.symMinBin, params.SYMBOL.symMaxBin, params.SYMBOL.symSeqLength, params.SYMBOL.symStdNum);
        end    
end

if isempty(k), 
    if isfield(params,'NumFeat') && ~isempty(params.NumFeat), 
        k = params.NumFeat; 
    else
        k = size(Y,2);
    end
else
    if k==-1, k=size(Y,2); end
end
if VERBOSE; fprintf(' Selecting %g features using %s algorithm (FEAST toolbox) ... ',k, params.MethodStr); end
tic;
switch params.Method

    case {1,2,3,4,5,6,7,8,9}
        ind = feast(params.MethodStr, k ,Y,label);
    case 10
        if ~iscell(params.MethodParams)
            MethodParams{1} = params.MethodParams;
        else
            MethodParams = params.MethodParams;
        end
        NumBeta = numel(MethodParams{1});
        ind = zeros(k,NumBeta);
        for i = 1:NumBeta
            ind(:,i) = feast( params.MethodStr, k, Y, label, MethodParams{1}(i) );
        end
    case 11
        P = allcomb(params.MethodParams);
        NumParams = size(P,1);
        ind = zeros(k,NumParams);
        for i = 1:NumParams
            ind(:,i) = feast( params.MethodStr, k, Y, label, P(i,1), P(i,2) );
        end
    case 12
        
end

if nargout == 2
    % Create weight vector
    W = (k:-1:1)'./k;
    R = zeros(n,size(ind,2));

    for i=1:size(ind,2)
        R(ind(:,i)) = W;
    end
else
    R=[];
end

end