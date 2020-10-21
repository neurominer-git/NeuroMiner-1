function [L, targscale, Lmin, Lmax, Lfact, Lpolyfact, Llogar] = nk_LabelTransform(Param, MODEFL, L)

targscale = 0; Lmin = min(L); Lmax = max(L); Lfact = 1; Lpolyfact = []; Llogar = [];

if isfield(Param,'LABELMOD') && strcmp(MODEFL,'regression')
    if isfield(Param.LABELMOD,'TARGETSCALE') && Param.LABELMOD.TARGETSCALE 
        fprintf('\n* Scaling of target labels [0<->1].')
        [L, IN] = nk_PerfScaleObj(L); L=full(L); Lmin = IN.minY; Lmax = IN.maxY; 
        targscale = 1;
    end
    if isfield(Param.LABELMOD,'POLYNOM') && ~isempty(Param.LABELMOD.POLYNOM)
        fprintf('\n* Computing polynomial transformation of target labels [labels.^%g].', Param.LABELMOD.POLYNOM)
        %Lmin = min(L); Lmax = max(L);
        L = L .^ Param.LABELMOD.POLYNOM;
        Lpolyfact = Param.LABELMOD.POLYNOM;
    end
    if isfield(Param.LABELMOD,'LOGAR') && ~isempty(Param.LABELMOD.LOGAR) && Param.LABELMOD.LOGAR
        fprintf('\n* Computing logarithm of target labels.')
        L = log(L); Llogar=1;
    end
end

end