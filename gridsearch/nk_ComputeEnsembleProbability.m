% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FORMAT: results = nk_ComputeEnsembleProbability(predictions, label)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Nikolaos Koutsouleris, 4/2012
function results = nk_ComputeEnsembleProbability(predictions, label, noscale, optcutoff, optcutoffperc)

global SCALE MODEFL RAND EVALFUNC MULTILABEL

if isempty(MODEFL), MODEFL = nk_input('Select prediction framework', 0, 'classification|regression'); end
if isempty(MULTILABEL), MULTILABEL.dim = size(predictions,2); end
if ~exist('label','var'), label = []; end
if ~exist('noscale','var'), noscale=[]; end
if ~exist('optcutoff','var'), optcutoff=[]; end
if ~exist('optcutoffperc','var'), optcutoffperc=[]; end
if isempty(EVALFUNC)
     switch MODEFL

        case 'regression'
            GridParam = nk_input('Cost function',0,'m', ...
                ['Mean Squared Error (MSE)|' ...
                'Normalized Root of Mean Squared Deviation (NRMSD)|' ...
                'Squared Correlation Coefficient (R^2)|' ...
                'Correlation coefficient'],[9 11 10 16],1);

        case 'classification'
            GridParam = nk_input('Cost function',0,'m', ...
                ['Accuracy|' ...
                'False Positive Rate|' ...
                'Positive Predictive Value|' ...
                'Matthews Correlation Coefficient|' ...
                'AUC|' ...
                'Gmean|' ...
                'Balanced Accuracy: BAC = ( sensitivity + specificity ) / 2)|' ...
                'Enhanced Balanced Accuracy: EBAC = ( sensitivity * specificity ) / 100)|' ...
                'F-score',],[1 4 5 6 7 13 14 17 15],1);

     end
     EVALFUNC = nk_DefineEvalFunc(GridParam); 
     
end

lx      = length(predictions);
results.num_predictions = zeros(lx,MULTILABEL.dim);

% Index to ~nan predictions
results.index_predictions = false(lx,MULTILABEL.dim);

% Mean prediction
results.mean_predictions = nan(lx,MULTILABEL.dim);

% STD of prediction
results.std_predictions = zeros(lx,MULTILABEL.dim);

% Confidence intervals of predictions
results.CI1_predictions = zeros(lx,MULTILABEL.dim); 
results.CI2_predictions = zeros(lx,MULTILABEL.dim);

if ~isempty(label)
    switch MODEFL
        case 'regression'
            results.R2          = zeros(1,MULTILABEL.dim);
            results.r           = zeros(1,MULTILABEL.dim);
            results.p           = zeros(1,MULTILABEL.dim);
            results.r_95CI_low  = zeros(1,MULTILABEL.dim);
            results.r_95CI_up   = zeros(1,MULTILABEL.dim);
            results.t           = zeros(1,MULTILABEL.dim);
        case 'classification'
    end
end
if isempty(noscale)
    [~, targscale, minX, maxX, ~, polyfact, logar] = nk_LabelTransform(SCALE, MODEFL, label);
else
    targscale=0; polyfact = []; logar = [];
end
if targscale, IN.revertflag = true; end

for curlabel=1:MULTILABEL.dim
    
    pred    = nan(lx,1); stdpred = pred; ci1 = pred; ci2 = pred;
    numpred = zeros(lx,1);
    mean_optcutoffs = nan(lx,1);
    std_optcutoffs = nan(lx,1);
    mean_optcutoffpercs = nan(lx,1);
    std_optcutoffpercs = nan(lx,1);
    %% Compute mean predictions and descriptive stats across CV2 permutations
    for i=1:lx % Loop through subjects

        predi = predictions{i,1,curlabel};
        if ~isempty(optcutoff)
            optcutoffi = optcutoff{i,1, curlabel};
            optcutoffperci = optcutoffperc{i,1, curlabel};
        end
            
        if ~isempty(predi)
            indnan = isnan(predi);
            if any(indnan)
                warning('NaN found in ensemble prediction. Some base model produced weird output!')
                if sum(indnan)==numel(predi), continue, end
            end
            predi           = predi(~indnan);
            if ~isempty(logar)
                predi       = exp(predi);
            end
            if ~isempty(polyfact), 
                predi       = predi.^(1/polyfact); end
            if targscale, 
                IN.minY = minX(curlabel); IN.maxY = maxX(curlabel);
                
                predi       = nk_PerfScaleObj(predi, IN); 
            end 

            lpredi          = numel(predi);
            numpred(i)      = lpredi;
            pred(i)         = nm_nanmedian(predi);
            stdpred(i)      = nm_nanstd(predi);
            try
                ci          = percentile(predi(~indnan),[2.5 97.5]);
            catch
                ci          = [1 1];
            end
            ci1(i)          = ci(1); 
            ci2(i)          = ci(2);
            if ~isempty(optcutoff), 
                mean_optcutoffs(i) = nm_nanmedian(optcutoffi); 
                std_optcutoffs(i) = nm_nanstd(optcutoffi); 
                mean_optcutoffpercs(i) = nm_nanmedian(optcutoffperci); 
                std_optcutoffpercs(i) = nm_nanstd(optcutoffperci); 
            end
        end
    end

    if isempty(label)
        ind = ~isnan(pred);
    else
        ind = ~isnan(pred) & ~isnan(label(:,curlabel));
        xlabl = label(ind,curlabel);
    end

    % Number of predictors per each hold-out sample
    results.num_predictions(:,curlabel) = numpred;

    % Index to ~nan predictions
    results.index_predictions(:,curlabel) = ind;

    % Mean prediction
    results.mean_predictions(:,curlabel) = pred;

    % STD of prediction
    results.std_predictions(:,curlabel) = stdpred;

    % Confidence intervals of predictions
    results.CI1_predictions(:,curlabel) = ci1; 
    results.CI2_predictions(:,curlabel) = ci2;
    
    % Cutoff for binarizing predicted risks of Cox-PH models
    if ~isempty(optcutoff),
        results.mean_cutoff_probabilities(:,curlabel) = mean_optcutoffs;
        results.std_cutoff_probabilities(:,curlabel) = std_optcutoffs;
        results.mean_cutoff_percentiles(:,curlabel) = mean_optcutoffpercs;
        results.std_cutoff_percentiles(:,curlabel) = std_optcutoffpercs;
    end
    
    switch MODEFL

        case 'regression'

             % Add grouping vector, if available
            if isfield(RAND,'groups'), results.grouping = RAND.groups; end

            %% Compute regression model statistics
            if ~isempty(label)
                results.costfun_crit(curlabel) = feval(EVALFUNC, xlabl, pred(ind));
                if numel(pred(ind))>1,
                    [r, p, rl, ru]                  = corrcoef(xlabl,pred(ind));
                    results.R2(curlabel)            = SCC(xlabl,pred(ind));
                    results.r(curlabel)             = r(1,2);
                    results.p(curlabel)             = p(1,2);
                    results.r_95CI_low(curlabel)    = rl(1,2);
                    results.r_95CI_up(curlabel)     = ru(1,2);
                    %% Convert correlation coefficient to T value
                    results.t(curlabel)             = r(1,2) * sqrt((sum(ind)-2) / (1 - r(1,2)*r(1,2)));
                else
                    results.R2(curlabel)            = NaN;
                    results.r(curlabel)             = NaN;
                    results.p(curlabel)             = NaN;
                    results.r_95CI_low(curlabel)    = NaN;
                    results.r_95CI_up(curlabel)     = NaN;
                    results.t(curlabel)             = NaN;
                end
                results.MAE(curlabel)               = MAERR(xlabl,pred(ind));
                results.MSE(:,curlabel)             = MSE(xlabl,pred(ind));
                results.NRSMD(:,curlabel)           = NRMSD(xlabl,pred(ind));
            end

        case 'classification'

            % These are the OOT predictions based on the average of decision values
            % / target predictions of the classifier ensembles

            %% Estimate OOT-probabilities
            % We would also like to compute the class probabilities of the 
            % ensembles majority vote as this is a more straight-forward way to
            % estimate the stability of OOT-classification for each test subject
            % only for models that do not produce risk estimates such
            % Cox-PH models
            if isempty(optcutoff),
                ul      = [1, -1];
                mx      = 2;
                prob    = nan(lx,mx);
                probmax = nan(lx,1);

                for i=1:lx
                    if ~ind(i), continue, end
                    spredi = sign(predictions{i,curlabel});
                    lpredi = length(spredi);
                    for j=1:mx
                        prob(i,j) = sum(spredi==ul(j))/lpredi;
                    end
                    [dum,ix] = max(prob(i,:));
                    lu = sum(prob(i,:) == dum);
                    if lu > 1 % throw coin if prediction is 50-50
                        lr = randperm(2,1);probmax(i) = ul(lr);
                    else
                        probmax(i) = ul(ix);
                    end

                end
                xprobmax = probmax(ind);

                %% Compute confusion matrix and classification performances based on OOT-probabilities
                results.prob_predictions = prob;
                results.prob_finalpred = xprobmax;
            end
            
            if ~isempty(label)
                if ~isempty(optcutoff),
                    I1 = pred > mean_optcutoffs; I2 = pred <= mean_optcutoffs;
                    results.prob_predictions = [pred 1-pred];
                    results.prob_finalpred = nan(lx,1);
                    results.prob_finalpred(I1)=1; results.prob_finalpred(I2)=0;
                    pred(ind) = pred(ind) - mean_optcutoffs(ind);
                    %% Compute confusion matrix based on optimally centered Cox-PH probabilities
                    results.costfun_crit = feval(EVALFUNC, xlabl, pred(ind));
                    results.contigency = ALLPARAM(xlabl, pred(ind)); 
                    results.prob_contigency = results.contigency;
                    results.globalcutoff_probabilities = mean_optcutoffs;
                    results.mean_globalcutoff_probabilities = nm_nanmean(mean_optcutoffs);
                    results.std_globalcutoff_probabilities = nm_nanstd(mean_optcutoffs);
                    results.mean_globalcutoff_percentiles = nm_nanmean(mean_optcutoffpercs);
                    results.std_globalcutoff_percentiles = nm_nanstd(mean_optcutoffpercs);
                 
                else
                    %% Compute confusion matrix based on decision values 
                    results.costfun_crit = feval(EVALFUNC, xlabl, pred(ind));
                    results.contigency = feval('ALLPARAM', xlabl, pred(ind)); 
                    results.prob_contigency = feval('ALLPARAM', xlabl, prob(ind,1)-0.5);  
                end
            end
            
    end
end

end