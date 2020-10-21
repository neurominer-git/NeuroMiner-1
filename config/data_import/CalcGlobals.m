function IO = CalcGlobals(IO)

%-Compute as mean voxel value (within per image fullmean/8 mask)
fprintf('Calculating globals\n')
IO.g = zeros(numel(IO.VV),1);
for i=1:numel(IO.VV), IO.g(i) = spm_global(IO.V(i)); fprintf('.'); end
if ~isfield(IO,'globnorm') || isempty(IO.globnorm), IO.globnorm = 1; end
