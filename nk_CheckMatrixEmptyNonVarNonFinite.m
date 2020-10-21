function [SIG, ind] = nk_CheckMatrixEmptyNonVarNonFinite(Y)

indany      =       any(Y);
indvar      =       std(Y) > 0;
indfinite   =       any(isfinite(Y));
ind         =       ~indany & ~indvar & ~indfinite;
if any(ind) , 
    warning('No-variance / empty / non-finite features in matrix'); 
    SIG = 1;
else
    SIG = 0;
end

end