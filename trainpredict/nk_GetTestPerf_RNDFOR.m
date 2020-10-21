function [rs, ds] = nk_GetTestPerf_RNDFOR(~, tXtest, ~, md, ~, ~)
global MODEFL
    
    switch MODEFL
        case 'classification'
            [rs, votes] = classRF_predict(tXtest,md); 
            ds = votes(:,2)./sum(votes,2);
        case 'regression'
            rs = regRF_predict(tXtest,md); ds=rs;
    end
end
