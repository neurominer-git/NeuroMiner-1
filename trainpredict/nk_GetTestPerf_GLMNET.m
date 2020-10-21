function [rs, ds] = nk_GetTestPerf_GLMNET(~, tXtest, ~, md, ~, ~)
global MODEFL

try
    p = glmnetPredict(md,tXtest);
catch
    tXtest
end
switch MODEFL
    case 'regression'
        rs = p(:,end); ds = p(:,end);
    case 'classification'
        rs = sign(p(:,end)-0.5) ; ds = p(:,end)-0.5;
end

