function [rs, ds] = nk_GetTestPerf_matLRN(~, tXtest, ~, md, ~, ~)
global MODEFL

p = md.predict(md, tXtest);
switch MODEFL
    case 'regression'
        rs = p; ds = p;
    case 'classification'
        if isstruct(p)
            rs = p.yhat; 
            if isfield(p,'D')
                ds = p.D;
            elseif isfield(p,'prob') 
                ds = p.prob-0.5;
            else
                ds = p.yhat;
            end
        else
            rs = p ; ds = p;
        end
end

