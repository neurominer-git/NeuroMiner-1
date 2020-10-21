function [DScore, idxF] = gflipmain(Y, labels, cpus, FeatInd, SortInd)

global FEATSEL
ix = size(Y,2);
if nargin < 6, if ~isempty(SortInd), Y = resamp(Y, SortInd, cpus);  end; end
if nargin < 4 || isempty(FeatInd), FeatInd = true(ix,1); end

gpu = FEATSEL.gflip.gpu; extra_param = FEATSEL.gflip.extra_param;

idxF    = false(ix,1);
FI      = find(FeatInd);
tY      = Y(:,FI);
switch gpu
    case 0
       [idx, DScore] = gflip(tY, labels, extra_param);
    case 1
       [idx, DScore] = gflip_gpu(tY, labels, extra_param);
end
idxF(FI(idx)) = true;

return
