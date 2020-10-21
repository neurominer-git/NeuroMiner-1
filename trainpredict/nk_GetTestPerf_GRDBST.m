function [rs, ds] = nk_GetTestPerf_GRDBST(~, tXtest, ~, md, ~, ~)
global MODEFL

ds = SQBMatrixPredict( md, single(tXtest));

switch MODEFL
    case 'regression'
        rs = ds; 
    case 'classification'
        rs = sign(ds);
end

