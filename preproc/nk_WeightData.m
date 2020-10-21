function [ wY, wYnew, wYoocv ] = nk_WeightData(W, TRd, CVd, OOCVd, SelectFlag, ScaleFlag)

wYnew = []; wYoocv = [];
if ~exist('ScaleFlag','var'), ScaleFlag = false; end
if size(W,2) == 1, W = W'; end
switch SelectFlag
    case 1
        ind = W~=0;
        fprintf(' Selecting %g non-zero voxels ...', sum(ind))
        wY = TRd.Y(:,ind);
        if exist('CVd','var') && ~isempty(CVd), wYnew = CVd.Y(:,ind); end
        if exist('OOCVd','var') && ~isempty(OOCVd),  wYoocv = OOCVd.Y(:,ind); end
    case 2
        fprintf(' Weighting data ...')
        wY = TRd.Y .* repmat(W, size(TRd.Y,1), 1);
        if exist('CVd','var') && ~isempty(CVd)
            wYnew = CVd.Y .* repmat(W, size(CVd.Y,1), 1);
        end
        if exist('OOCVd','var') && ~isempty(OOCVd)
            wYoocv = OOCVd.Y.* repmat(W, size(OOCVd.Y,1), 1);
        end
end

if ScaleFlag
    IN.AcMatFl = true;
    [wY, IN] = nk_PerfScaleObj(wY,IN);
    if exist('CVd','var') && ~isempty(CVd)
        wYnew = nk_PerfScaleObj(wYnew, IN);
    end
    if exist('OOCVd','var') && ~isempty(OOCVd)
        wYoocv = nk_PerfScaleObj(wYoocv, IN);
    end
end

end