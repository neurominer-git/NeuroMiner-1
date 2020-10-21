function r = rfe_algo_settings(Y, label, Ynew, labelnew, Ps, FullFeat, FullParam, ActStr)

global RFE VERBOSE SVM TRAINFUNC

r.FullInd = find(any(Y) & std(Y)~=0 & sum(isfinite(Y))~=0 & FullFeat==1); r.Y = Y(:,r.FullInd); r.Ynew = Ynew(:,r.FullInd);
r.kFea = numel(r.FullInd); 
r.FeatRandPerc = 0;
r.FeatStepPerc = true;

switch RFE.Wrapper.type
    case 1
        %% Feature Sorting Mode
        % if WeightSort==2 then the feature weight vector will be used for sorting
        % the features. This is only valid for linear models
        if ~isfield(RFE.Wrapper.GreedySearch,'WeightSort') || isempty(RFE.Wrapper.GreedySearch.WeightSort)
            r.WeightSort = 1; 
        else
            r.WeightSort = RFE.Wrapper.GreedySearch.WeightSort;
        end

        %% Feature block size settings
        kFea = length(r.FullInd); 
        if ~isfield(RFE.Wrapper.GreedySearch,'FeatStepPerc')
            r.lperc = 1/(kFea/100); % lperc should multiply to a feature stepping of 1 
            r.FeatStepPerc = false;
        elseif ~RFE.Wrapper.GreedySearch.FeatStepPerc
            r.lperc = 1;
        else
            r.lperc = RFE.Wrapper.GreedySearch.FeatStepPerc;
        end

        %% Early stopping setting
        r.MinNum = 1;
        if RFE.Wrapper.GreedySearch.EarlyStop.Thresh
           switch RFE.Wrapper.GreedySearch.EarlyStop.Perc
               case 1 % Percentage Mode
                   r.MinNum = size(r.Y,2) / 100 * RFE.Wrapper.GreedySearch.EarlyStop.Thresh; 
               case 2 % Absolute mMde
                   r.MinNum = RFE.Wrapper.GreedySearch.EarlyStop.Thresh;
           end
        end
        if isfield(RFE.Wrapper.GreedySearch,'FeatRandPerc')
            r.FeatRandPerc = RFE.Wrapper.GreedySearch.FeatRandPerc;
        end
end

% Optimization criterion setup
[r.evaldir, ~, r.optfunc, r.optparam] = nk_ReturnEvalOperator(SVM.GridParam);

%% Obtain full feature space performance, if not already done so
if ~exist('FullParam','var') || isempty(FullParam)
    
    switch ActStr
        case 'Tr'
            r.T = r.Y;
            r.L = label;
        case 'CV'
            r.T = r.Ynew;
            r.L = labelnew;
        case 'TrCV'
            r.T = [r.Y; r.Ynew];
            r.L = [label; labelnew];
    end
    
    [~, r.FullModel] = feval(TRAINFUNC, r.Y, label, 1, Ps);    
    r.FullParam = nk_GetTestPerf(r.T, r.L, [], r.FullModel, r.Y);
    if VERBOSE, fprintf('\nFull model:\t# Features: %g, %s = %g', ...
                        numel(r.FullInd), ActStr, r.FullParam), end
else
    r.FullParam = FullParam;
    r.FullModel = [];
end

if isfield(RFE.Wrapper.GreedySearch,'KneePointDetection') && RFE.Wrapper.GreedySearch.KneePointDetection == 1
    r.KneePoint = true;
else
    r.KneePoint = false;
end

if isfield(RFE.Wrapper.GreedySearch,'CaseFrwd'),
    r.CaseU = RFE.Wrapper.GreedySearch.CaseFrwd.Upper;
    r.CaseL = RFE.Wrapper.GreedySearch.CaseFrwd.Lower;
    r.CaseStep = RFE.Wrapper.GreedySearch.CaseFrwd.Step;
    r.CaseSpU = r.CaseU(2):-1*round((r.CaseU(2)-r.CaseU(1))/r.CaseStep):r.CaseU(1);
    r.CaseSpL = r.CaseL(1):round((r.CaseL(2)-r.CaseL(1))/r.CaseStep):r.CaseL(2);
end
