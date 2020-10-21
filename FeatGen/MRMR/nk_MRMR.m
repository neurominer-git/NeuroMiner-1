function ind = nk_MRMR(Y, label, k, params)

global MODEFL VERBOSE

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

if isempty(k), k = params.NumFeat; end

switch params.mRMRscheme
    case 1
        ind = mrmr_mid_d(Y,label,k); % This is the MID scheme
    case 2
        ind = mrmr_miq_d(Y,label,k); % This is the MIQ scheme
end

end