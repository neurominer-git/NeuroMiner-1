function I = nk_VisModels_ExtraLabels(I, inp, nperms, compfun)
global EVALFUNC

nL = size(inp.extraL.L,2);
I.EXTRA_L = struct('VCV2MORIG_EVALFUNC_GLOBAL', nan(inp.nclass,1), 'VCV2MPERM_EVALFUNC_GLOBAL', nan(inp.nclass,nperms(1)), 'VCV2MPERM_GLOBAL', nan(inp.nclass,nperms(1)));

for h=1:inp.nclass
    
    indempt = cellfun(@isempty,I.VCV2MORIG_S);
    Porig = cellfun(@nm_nanmedian,I.VCV2MORIG_S(~indempt(:,h),h)); 
    indnonnan = ~isnan(Porig); 
    
    if inp.targscale, IN.revertflag = true; IN.minY = inp.minLbCV; IN.maxY = inp.maxLbCV; Porig = nk_PerfScaleObj(Porig, IN); end
        
    for l = 1:nL
        
        indnonnanL = ~isnan(inp.extraL.L(:,l)) & indnonnan;
        LO = inp.extraL.L(indnonnanL,l);
        PO = Porig(indnonnanL); 
        I.EXTRA_L(l).LABEL_NAME = inp.extraL.Lnames{l};
        I.EXTRA_L(l).LABEL = inp.extraL.L(:,l);
        I.EXTRA_L(l).VCV2MORIG_EVALFUNC_GLOBAL = feval(EVALFUNC, LO, PO);
        I.EXTRA_L(l).CONTINGENCY = ALLPARAM(LO,PO);
        fprintf('\nExternal Label %s: Testing observed model performance against permuted models using entire data: %g permutations\n', I.EXTRA_L(l).LABEL_NAME, nperms(1))
        for perms = 1 : nperms(1)
            PP = cellfun(@nm_nanmedian, I.VCV2MPERM_S(~indempt(:,h),h,perms));
            if inp.targscale, PP = nk_PerfScaleObj(PP, IN); end
            I.EXTRA_L(l).VCV2MPERM_EVALFUNC_GLOBAL(h,perms) = feval(EVALFUNC, LO, PP(indnonnanL)); 
            crt = feval(compfun, I.EXTRA_L(l).VCV2MPERM_EVALFUNC_GLOBAL(h,perms), I.EXTRA_L(l).VCV2MORIG_EVALFUNC_GLOBAL);
            if ~crt, fprintf('*'); else, fprintf('.'); end
            I.EXTRA_L(l).VCV2MPERM_GLOBAL(h,perms) = crt;
        end
        I.EXTRA_L(l).P = sum(I.EXTRA_L(l).VCV2MPERM_GLOBAL(h,:))/nperms(1);
        fprintf('\nP value: %g',I.EXTRA_L(l).P );
    end

end