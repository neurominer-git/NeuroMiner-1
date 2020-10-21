function GDanalysis = nk_MultiPerfComp(GDanalysis, multi_pred, label, ngroups, act)

if ~exist('act','var') || isempty(act), act = 'pred'; end  
switch act
    case 'pred'
        fld = 'MultiClass'; 
        % Compute multi-class probabilties
        multi_prob = nk_ConvProbabilities(multi_pred, ngroups);
    case 'prob'
        fld = 'MultiClassProb'; 
        multi_prob = multi_pred;
end

lx = size(multi_prob,1);
pred = nan(lx,1); 
stdpred = pred; ci1 = pred; ci2 = pred;

% Loop through cases
if iscell(multi_pred)
    for i=1:lx
        if isempty(multi_pred{i})
            pred(i) = NaN; stdpred(i) = NaN; ci1(i) = NaN; ci2(i) = NaN; errs(i) = NaN;
        else
            % Maximum probability decides about multi-class membership
            [maximum,pred(i)] = max(multi_prob(i,:),[],'includenan'); 
            if isnan(maximum), pred(i)=NaN; continue; end
            % Is this useful: ?
            stdpred(i) = std(multi_pred{i});
            ci = percentile(multi_pred{i},[2.5 97.5]);
            ci1(i) = ci(1); ci2(i) = ci(2);
        end
    end
else
    %numpred = size(multi_pred,2);
    [maximum,pred] = max(multi_prob,[],2,'includenan'); 
    inan = isnan(maximum);
    pred(inan)=NaN;
    stdpred = std(multi_pred,[],2);
    stdpred(inan)=NaN;
    ci = cell2mat(arrayfun( @(i) percentile(multi_pred(i,:),[2.5 97.5]),1:lx,'UniformOutput',false )');
    ci(inan)=NaN;
    ci1 = ci(:,1); ci2 = ci(:,2);
end

% Compute multi-class performance
GDanalysis = nk_MultiEvalPerf(GDanalysis, label, pred, ngroups, fld);

% Store prediction results in structure
GDanalysis.(fld).multi_probabilitiesCV2     = multi_prob;
GDanalysis.(fld).multi_predictionsCV2       = pred;
GDanalysis.(fld).multi_predictionsCV2_std   = stdpred;
GDanalysis.(fld).multi_predictionsCV2_ci1   = ci1;
GDanalysis.(fld).multi_predictionsCV2_ci2   = ci2;

end