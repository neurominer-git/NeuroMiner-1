function inp = nk_PerfInitSpatial(analysis, inp, paramfl)
global VERBOSE MODEFL

kbin = inp.nclass;
nM  = numel(inp.PREPROC);
Yocv = []; Cocv = [];

for n=1:nM
    
    if nM>1
        Y = inp.X(n).Y; 
        if isfield(inp.X(n),'Yocv') && ~isempty(inp.X(n).Yocv), Yocv = inp.X(n).Yocv; end
        PREPROC = inp.PREPROC{n};
    else
        Y = inp.X.Y; 
        if isfield(inp.X(n),'Yocv') && ~isempty(inp.X(n).Yocv), Yocv = inp.X.Yocv; end 
        PREPROC = inp.PREPROC; 
    end
    
    if iscell(paramfl)
        tparamfl = paramfl{n};
    else
        tparamfl = paramfl;
    end
    
    tparamfl.PV = inp.X(n);
    if VERBOSE, fprintf('\nGenerate pre-processing parameter array for CV2 partition [%g,%g].\n',inp.f,inp.d); end
    tparamfl = nk_PrepPreprocParams(PREPROC, tparamfl, analysis, n, inp.ll, inp.curlabel);
    
    if iscell(PREPROC)
        BINMOD = PREPROC{1}.BINMOD;
    else
        BINMOD = PREPROC.BINMOD;
    end
    if BINMOD || strcmp(MODEFL,'regression')
        ukbin = kbin;   if VERBOSE; fprintf('\nProcessing Mode: binary / regression preprocessing'); end
    else
        ukbin = 1;      
        if VERBOSE; fprintf('\nProcessing Mode: multi-group preprocessing'); end
    end
    
    [ sY, sYocv, sCocv, inp ] = nk_PerfPreprocessSpatial( Y, Yocv, Cocv, inp, tparamfl, BINMOD, kbin, ukbin);
    
    if nM>1 
        inp.X{n}.sY = sY;
        if ~isempty(sYocv), inp.X{n}.sYocv = sYocv; end
        if ~isempty(sCocv), inp.X{n}.sCocv = sCocv; end
    else
        inp.X(n).sY = sY;
        if ~isempty(sYocv), inp.X(n).sYocv = sYocv; end
        if ~isempty(sCocv), inp.X(n).sCocv = sCocv; end
    end
    

end