function [Thresh, Crit, Mx, YnewPred, YPerf, YnewPerf] = nk_DefineThresh(FEATSEL, Y, DScore, labels, Ynew, labelnew, sigma, curclass)

Mx = []; YnewPred = []; YPerf=[]; YnewPerf=[];
% Determine threshold 
switch FEATSEL.salthreshmode
    case 0
        Thresh = []; Crit = [];
    case 1 % Percentile
        if numel(FEATSEL.salCI) > 1
            Thresh = percentile(DScore, FEATSEL.salCI(curclass));
            Crit = FEATSEL.salCI(curclass);
        else
            Thresh = percentile(DScore, FEATSEL.salCI);
            Crit = FEATSEL.salCI;
        end
    case 2 % Absolute
        if numel(FEATSEL.salThresh) > 1
            Thresh = FEATSEL.salThresh(curclass);
            Crit = FEATSEL.salThresh(curclass);
        else
            Thresh = FEATSEL.salThresh;
            Crit = FEATSEL.salThresh;
        end
        
    case 3 % Adaptive
        [Thresh, Crit, Mx, YnewPred, YPerf, YnewPerf] = nk_AdaptiveThreshold2(Y, labels, DScore, FEATSEL.salAdap, Ynew, labelnew, sigma);
end

end