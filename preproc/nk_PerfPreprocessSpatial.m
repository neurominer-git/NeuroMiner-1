function [sY, sYocv, sCocv, inp, optfl, ukbin, uBINMOD, BINMOD] = ...
    nk_PerfPreprocessSpatial( Y, Yocv, Cocv, inp, paramfl, BINMOD, kbin, ukbin)

global PREPROC
sYocv = []; sCocv = []; 

%Eventually, apply spatial operations on image data
if isfield(paramfl,'PREPROC') && isfield(paramfl.PREPROC,'SPATIAL') && paramfl.PREPROC.SPATIAL.cubetype>1 
    % Filter the data (imaging only)
    smoothfl = true; if isfield(inp,'issmoothed') && inp.issmoothed, smoothfl = false; end
    sY = cell(kbin,1); 
    if ~isempty(Yocv), sYocv = cell(kbin,1); end
    if ~isempty(Cocv), sCocv = cell(kbin,1); end
    optfl = true; ukbin = kbin; 
    uBINMOD = 0; 
    if smoothfl
        if inp.multiflag
            uPREPROC = nk_SetParamChain(paramfl, 1, PREPROC);
            tP = nk_ReturnParamChain(uPREPROC, 1); 
            fprintf('\nAll Predictors: Smoothing Training & CV data');
            tsY = nk_PerfSpatFilt2( Y, uPREPROC, paramfl.PV ); 
            if isfield(inp.X,'Yw'), inp.Yw = cell(kbin,1); tsYw = nk_PerfSpatFilt2( inp.X.Yw, uPREPROC, paramfl.PV ); end
            if ~isempty(Yocv), tsYocv = nk_PerfSpatFilt2( Yocv, uPREPROC, paramfl.PV ); end 
            if ~isempty(Cocv), tsCocv = nk_PerfSpatFilt2( Cocv, uPREPROC, paramfl.PV ); end
            for u=1:kbin
                paramfl.P{u} = tP;
                sY{u} = tsY;
                if isfield(inp.X,'Yw'), inp.Yw{u} = tsYw; end
                if ~isempty(Yocv), sYocv{u} = tsYocv; end
                if ~isempty(Cocv), sCocv{u} = tsCocv; end
            end
        else
            for u=1:kbin
                uPREPROC = nk_SetParamChain(paramfl, u, PREPROC);
                paramfl.P{u} = nk_ReturnParamChain(uPREPROC, 1); 
                fprintf('\nPredictor #%g: Smoothing Training & CV data',u)
                sY{u} = nk_PerfSpatFilt2( Y, uPREPROC, paramfl.PV ); 
                if isfield(inp.X,'Yw'), 
                    fprintf('\nSmoothing weighting map')
                    inp.Yw{u} = nk_PerfSpatFilt2( inp.X.Yw, uPREPROC, paramfl.PV ); 
                else
                    I = arrayfun( @(j) isfield(uPREPROC.ACTPARAM{j},'RANK'), 1:numel( uPREPROC.ACTPARAM ));
                    if any(I), 
                        inp.Yw{u} = uPREPROC.ACTPARAM{I}.RANK.EXTERN;
                        inp.Yw{u} = nk_PerfSpatFilt2( inp.Yw{u}, uPREPROC, paramfl.PV ); 
                    end
                end
                
                
                if ~isempty(Yocv), 
                    fprintf('\nPredictor #%g: Smoothing independent test data (%s)',u, inp.desc_oocv);
                    sYocv{u} = nk_PerfSpatFilt2( Yocv, uPREPROC, paramfl.PV ); 
                end
                if ~isempty(Cocv), 
                    fprintf('\nPredictor #%g: Smoothing calibration data',u)
                    sCocv{u} = nk_PerfSpatFilt2( Cocv, uPREPROC, paramfl.PV ); 
                end
            end
        end
    else
        uBINMOD = BINMOD;
        sY = Y;
        if ~isempty(Yocv), sYocv = Yocv; end
        if ~isempty(Cocv), sCocv = Cocv; end
        if isfield(inp,'X') && isfield(inp.X,'Yw') && ~isfield(inp,'Yw'), inp.Yw = inp.X.Yw; end
    end
else
    uBINMOD = BINMOD;
    optfl = false;
    sY = Y;
    if ~isempty(Yocv), sYocv = Yocv; end
    if ~isempty(Cocv), sCocv = Cocv; end
    if isfield(inp,'X') && isfield(inp.X,'Yw') && ~isfield(inp,'Yw'), inp.Yw = inp.X.Yw; end
end