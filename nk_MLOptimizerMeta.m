function [ METAres, IN ] = nk_MLOptimizerMeta( analysis, IN )

global FUSION META MODEFL 

GDdims = analysis.GDdims;
aggmode = 1; METAres = [];
% Perform 2nd-level learning
if FUSION.flag == 3 && ~isempty(META)
    switch META.flag 
        case 1 % Bagging
            ind0 = find(~cellfun(@isempty,GDdims));
            nG = numel(ind0);
            P = cell(IN.nclass,1);
            if isfield(GDdims{ind0(1)},'multi_predictions')
               mP = cell(size(IN.labels,1),1);
               best_MultiCVperf = zeros(nG,1);
               best_MultiTSperf = zeros(nG,1);
            end
            switch MODEFL
                case 'classification'
                    Perf            = zeros(nG,IN.nclass);
                case 'regression'
                    R2      = zeros(nG,1);
                    R       = zeros(nG,1);
                    T       = zeros(nG,1);
                    NRSMDs  = zeros(nG,1);
                    Perf    = zeros(nG,1);
            end
            best_CVperf     = zeros(nG,IN.nclass);
            best_TSperf     = zeros(nG,IN.nclass);
            best_Complexity = zeros(nG,IN.nclass);
            best_Error      = zeros(nG,IN.nclass);
            
            switch aggmode
                case 1
                    % Here we aggregate only the median predictions for
                    % each case. That is the bagged predictor
                    for curclass=1:IN.nclass
                        P{curclass} = zeros(numel(IN.labels),nG);
                        for n=1:nG
                            nGDdims = GDdims{ind0(n)};
                            switch MODEFL
                                case 'classification'
                                    Pc                  = nGDdims.BinClass{curclass}.mean_predictions;
                                    Perf(n,curclass)    = nGDdims.BinClass{curclass}.costfun_crit;
                                    
                                case 'regression'
                                    Pc                  = nGDdims.Regr.mean_predictions;
                                    Perf(n)             = nGDdims.Regr.costfun_crit;
                                    R(n)                = nGDdims.Regr.r;
                                    R2(n)               = nGDdims.Regr.R2;
                                    T(n)                = nGDdims.Regr.t;
                                    NRSMDs(n)           = nGDdims.Regr.NRSMD;
                            end
                            P{curclass}(:,n)            = Pc;
                            if isfield(nGDdims,'multi_predictions')
                                if n==1
                                    mP                 = nGDdims.multi_predictions;
                                else
                                    mP                 = arrayfun( @(i) cellmat_mergecols(mP(i), nGDdims.multi_predictions(i)),1:size(IN.labels,1))';
                                end
                                best_MultiCVperf(n)         = nGDdims.best_MultiCVperf;
                                best_MultiTSperf(n)         = nGDdims.best_MultiTSperf;
                            end
                            best_CVperf(n,:)                = cell2mat(nGDdims.best_CVperf);
                            best_TSperf(n,:)                = cell2mat(nGDdims.best_TSperf);
                            best_Complexity(n,:)            = cell2mat(nGDdims.best_Complexity);
                            best_Error(n,:)                 = cell2mat(nGDdims.best_Error);
                        end
                        P{curclass}                     = mat2cell(P{curclass}, ones(size(P{curclass},1),1)); 
                    end
                case 2
                    % Here we aggregate all base learners' predictions
                    P = cell(IN.nclass,1); 
                    for curclass=1:IN.nclass
                        P{curclass}                     = GDdims{1}.predictions(:,curclass);
                        for n=2:nG
                            nGDdims                     = GDdims{ind0(n)};
                            P{curclass}                 = cellmat_mergecols(P{curclass}, nGDdims.predictions(:,curclass));
                        end
                    end
            end
    end     
    
    switch MODEFL
        case 'classification'
            lu = [1 -1];
            L = nan(numel(IN.labels),IN.nclass);
            for curclass=1:IN.nclass
                cvu = analysis.params.cv.class{1,1}{curclass};
                for y = 1:numel(cvu.groups)
                    indy = IN.labels == cvu.groups(y); L(indy,curclass) = lu(y);
                end
                if y==1, indy = IN.labels ~= cvu.groups(y);L(indy,curclass) = lu(y); end
                METAres.BinClass{curclass} = nk_ComputeEnsembleProbability(P{curclass}, L(:,curclass));
                best_TSperf_multimodal(curclass) = METAres.BinClass{curclass}.costfun_crit;
            end
            % Multi-group data available?
            if exist('mP','var')
                METAres = nk_MultiPerfComp(METAres, mP, IN.labels, IN.ngroups);
                METAres.best_MultiCVperf_unimodal = best_MultiCVperf;
                METAres.best_MultiTSperf_unimodal = best_MultiTSperf;
            end
        case 'regression'
            L = IN.labels; 
            METAres.Regr = nk_ComputeEnsembleProbability( P{1}, L , 1);
            METAres.r = R;
            METAres.t = T;
            METAres.R2 = R2;
            METAres.NRSMD = NRSMDs;
            METAres.best_TSperf_multimodal = METAres.Regr.costfun_crit;
    end
    METAres.labels = L;
    METAres.predictions = P;
    METAres.best_TSperf_multimodal = best_TSperf_multimodal;
    METAres.best_CVperf_unimodal = best_CVperf;
    METAres.best_TSperf_unimodal = best_TSperf;
    METAres.best_Complexity_unimodal = best_Complexity;
    METAres.best_Error_unimodal = best_Error;
end