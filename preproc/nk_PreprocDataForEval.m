function [ wY, wYnew, wYoocv ] = nk_PreprocDataForEval(param, wY, wYnew, wYoocv)

if isfield(param,'DR') && isfield(param.DR,'RedMode')
    fprintf('%s',param.DR.RedMode)
    [wY, mpp, indNonZero] = nk_PerfRed2(wY, param.DR);
    if exist('wYnew','var') && ~isempty(wYnew)
        wYnew = nk_PerfRed2(wYnew, param.DR, mpp, indNonZero);
    else
        wYnew = [];
    end
    if exist('wYoocv','var') && ~isempty(wYoocv)
        wYoocv = nk_PerfRed2(wYoocv, param.DR, mpp, indNonZero);
    else
        wYoocv = [];
    end
end

if isfield(param,'SCALE') && param.SCALE.flag
    [wY, IN] = nk_PerfScaleObj(wY);
    if exist('wYnew','var') && ~isempty(wYnew)
        wYnew = nk_PerfScaleObj(wYnew, IN);
    else
        wYnew = [];
    end
    if exist('wYoocv','var') && ~isempty(wYoocv)
        wYoocv = nk_PerfScaleObj(wYoocv, IN);
    else
        wYoocv = [];
    end
end

end

