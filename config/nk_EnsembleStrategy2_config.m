function param = nk_EnsembleStrategy2_config(param, SVM, MODEFL, defaultsfl, parentstr)

Ensemble.type           = 0;    % Aggregated ensemble
Ensemble.Metric         = 2;    % Decision values (if SVM), Probabilities if RVM
Ensemble.MinNum         = 1;    % Min # of classifiers
Ensemble.Perc           = 75;   % Min % of cross-subspace feature overlap
Ensemble.DataType       = 2;    % optimize according to CV1 test data prediction performance
Ensemble.Weighting      = 0;    % Weight base hypotheses of ensemble classifier
Ensemble.DivCrit        = 2;
Ensemble.ConstructMode  = 0;
Ensemble.DivStr         = '';
Ensemble.DivFunc        = '';
Ensemble.OptFunc        = '';
Ensemble.EntropRegMode  = 1;
SubSpaceStrategy        = 1;   % Winner takes all
SubSpaceCrit            = 0;   % No subspace threshold
act                     = 0;
CostFun                 = 2;

if ~isempty(SVM) && isfield(SVM,'GridParam')
    switch SVM.GridParam
        case {1,5,6,7,10,13,14,15,17}
            OPa = ' >= '; OPb = ' > ';
            Ensemble.CompFunc = 'max';
        otherwise
            OPa = ' <= '; OPb = ' < ';
            Ensemble.CompFunc = 'min';
    end
else
    warndlg('Setup of prediction algorithm and main optimization parameter is required!')
    return
end

%% Setup of ensemble construction strategy
if ~exist('defaultsfl','var') || isempty(defaultsfl),  defaultsfl = 0; end;

if ~defaultsfl 
    
    if isfield(param,'SubSpaceStrategy'), SubSpaceStrategy = param.SubSpaceStrategy; end
    if isfield(param,'SubSpaceCrit'), SubSpaceCrit = param.SubSpaceCrit; end
    if isfield(param,'CostFun'), CostFun = param.CostFun; end
    if isfield(param,'EnsembleStrategy') && isfield(param.EnsembleStrategy,'DataType'), 
        Ensemble = param.EnsembleStrategy; 
        if ~isfield(Ensemble,'EntropRegMode'), Ensemble.EntropRegMode = 1; end
    else
        param.EnsembleStrategy = Ensemble;
    end
    if ~CostFun, param.CostFun = 2; end
    
    d = nk_GetParamDescription2([],param,'EnsType'); menuact = [];
    
    menustr = sprintf('Select subspace selection strategy [ %s ]', d.SubSpaceStrat);                                                menuact = [ menuact 1 ];
    
    if SubSpaceStrategy > 1
        menustr = sprintf('%s|Define ensemble optimization method [ %s ]', menustr, d.EnsConMode);                                  menuact = [ menuact 2 ];
        if Ensemble.ConstructMode && Ensemble.ConstructMode ~= 4
            menustr = sprintf('%s|Choose optimization function [ %s ]', menustr, d.EnsDivCrit);                                     menuact = [ menuact 3 ];
        end
        if Ensemble.type ~= 0 && Ensemble.type ~= 9  
            if strcmp(MODEFL,'classification')
                menustr = sprintf('%s|Use algorithm output scores or label predictions [ %s ]', menustr, d.FeatMetric);             menuact = [ menuact 4 ];
            else
                Ensemble.Metric = 1; 
            end
            menustr = sprintf('%s|Define minimum # of classifier to be selected [ %1.0f ]', menustr, Ensemble.MinNum);              menuact = [ menuact 5 ];
            menustr = sprintf('%s|Enable weighting of feature subspaces [ %s ]', menustr, d.EnsWeighting);                          menuact = [ menuact 6 ];
        elseif Ensemble.type == 9  
            menustr = sprintf('%s|Define percentage of cross-subspace feature agreement [ %g ]', menustr, Ensemble.Perc);           menuact = [ menuact 7 ];
            menustr = sprintf('%s|Define minimum number of features to select across features subspaces [ %1.0f ]', ...
                                                                                                        menustr, Ensemble.MinNum);  menuact = [ menuact 8 ];
            Ensemble.Weighting = 0;
        end
    end
    
    nk_PrintLogo
    mestr = 'Subspace-based ensemble optimization setup'; navistr = [parentstr ' >>> ' mestr]; cprintf('*blue','\nYou are here: %s >>>',parentstr); 
    act = nk_input(mestr, 0,'mq', menustr, menuact);
    
    switch act

        case 1
             switch SVM.GridParam
                case {9,11,12,18}
                   mxstr ='minimum '; threshstr = 'BELOW'; df = 25;
                 otherwise
                   mxstr= 'maximum'; threshstr = 'ABOVE'; df = 75;
             end
             SubSpaceStrategy = nk_input('Subspace selection strategy',0, 'm', ...
                        ['Subspace with ' mxstr d.CostType ' criterion (winner takes it all)|' ...
                        'Subspace ensemble with ' mxstr d.CostType ' within a range of the maximum|' ...
                        'Subspace ensemble ' mxstr d.CostType ' ' threshstr ' a given percentile|' ...
                        'All-subspace ensemble'],1:4, SubSpaceStrategy);
             switch SubSpaceStrategy
                case 1
                    SubSpaceCrit = 0; Ensemble.ConstructMode = 0;
                case 2
                    SubSpaceCrit = nk_input(['Range from ' d.CostType mxstr],0,'e', 5);
                case 3
                    SubSpaceCrit = nk_input(['Percentile [%] for ' d.CostType ' cutoff'],0, 'e', df);
                case 4
                    SubSpaceCrit = 0; Ensemble.ConstructMode = 0;
             end
        case 2
             Ensemble.ConstructMode = nk_input('Select ensemble construction method',0,'m', ...
                   ['Simply aggregate all learners into ensemble|' ...
                    'Optimize ensemble using backward base learner elimination|' ...
                    'Optimize ensemble using forward base learner selection|' ...
                    'Create a single classifier using probabilistic feature subspace construction|' ...
                    'Use AdaBoost to determine optimal weighting of bases learners'],[0 1 2 4 5], Ensemble.ConstructMode);
             if (Ensemble.ConstructMode && Ensemble.ConstructMode ~= 4) && ~Ensemble.DivCrit, Ensemble.DivCrit = 2; end
        case 3
            Ensemble.DivCrit = nk_input('Select the optimization function',0,'m', ...
                    ['Optimize entropy of ensemble''s component models|' ...
                     'Optimize ensemble''s prediction performance & entropy|' ...
                     'Optimize Kappa-Diversity among ensemble''s component decisions|' ...
                     'Optimize Bias-Variance decomposition of ensemble''s component errors'],1:4, Ensemble.DivCrit);
             if Ensemble.DivCrit == 2
                 Ensemble.EntropRegMode = nk_input('Choose optimization function',0, 'm', ...
                                         ['Entropy Regularization: Prediction performance >= Param AND Entropy >= Param|' ...
                                          'Mixed criterion: (Prediction performance > Param AND Entropy >= Param) OR ' ...
                                          '(Prediction performance >= Param AND Entropy > Param)'],[1,2],2); 
             end
        case 4
            Ensemble.Metric = nk_input('Use predicted labels or algorithm scores for ensemble construction',0,'m', ...
                                        ['Predicted labels (Hard decision ensemble)|' ...
                                         'Algorithm scores (Soft decision ensemble)'], 1:2, Ensemble.Metric);
        case 5
            Ensemble.MinNum = nk_input('Minimum number of classifiers to retain', 0,'e',Ensemble.MinNum);
        case 6
            if ~Ensemble.Weighting, Ensemble.Weighting=2; end
            Ensemble.Weighting = nk_input('Weight base hypotheses?',0,'m', ...
                'No weighting|Weight = 1 / resubstitution error',[0,1], Ensemble.Weighting);
        case 7
            Ensemble.Perc = nk_input('[%] of cross-subspace feature agreement',0,'e',Ensemble.Perc);
        case 8
            Ensemble.MinNum = nk_input('Minimum number of features to retain',0,'e',Ensemble.MinNum);
    end
end

switch Ensemble.ConstructMode
    case 0
       Ensemble.type = 0;
    case 1
        if Ensemble.DivCrit == 5, Ensemble.type = 2; else Ensemble.type = Ensemble.DivCrit; end
    case 2
        if Ensemble.DivCrit == 5, Ensemble.type = 6; else Ensemble.type = Ensemble.DivCrit + 4; end
    case 4
        Ensemble.type = 9;
    case 5
        %% to be programmed
end

%% Define Ensemble-based subspace selection method
if SubSpaceStrategy > 1
    switch Ensemble.DivCrit
        case 0
            if strcmp(MODEFL,'classification')
                Ensemble.DivFunc = 'nk_Entropy';
            else
                Ensemble.DivFunc = 'nk_RegAmbig';
            end
        case {1,2}
            Ensemble.DivFunc = 'nk_Entropy'; 
            switch Ensemble.DivCrit
                case 1
                    Ensemble.DivStr = 'entropy'; Ensemble.OptFunc = 'EntropMax';
                case 2
                    Ensemble.DivStr = 'prediction performance (reg. entropy)';
                    Ensemble.OptFunc = 'CVMax';
                    % Define Inline Optimization function depending on
                    % Optimization Parameter and regularization mode
                    switch Ensemble.EntropRegMode
                        case 1
                            OptInline1 = ['Perf' OPa 'PerfCrit && Entropy >= EntropyCrit'];
                            OptInline2 = ['OrigPerf' OPb 'OptPerf || (OrigPerf == OptPerf && OrigEntropy > OptEntropy)'];
                        case 2
                            OptInline1 = ['(Perf' OPb 'PerfCrit && Entropy >= EntropyCrit) || (Perf' OPa 'PerfCrit && Entropy > EntropyCrit)'];
                            OptInline2 = ['OrigPerf' OPb 'OptPerf || OrigEntropy > OptEntropy'];
                    end
                    Ensemble.OptInlineFunc1  = inline(OptInline1, 'Perf', 'PerfCrit', 'Entropy', 'EntropyCrit');
                    Ensemble.OptInlineFunc2  = inline(OptInline2, 'OrigPerf', 'OptPerf', 'OrigEntropy', 'OptEntropy');
            end
        case 3
            Ensemble.DivFunc = 'nk_Diversity'; Ensemble.OptFunc = 'DiversityMax';
            Ensemble.DivStr = 'Kappa-Diversity';
        case 4
            Ensemble.DivFunc = 'nk_Lobag'; Ensemble.OptFunc = 'EDMin';
            Ensemble.DivStr = 'Bias-Variance decomposition of ensemble error';
        case 5 % for regression models
            Ensemble.DivFunc = 'nk_RegAmbig'; Ensemble.OptFunc = 'CVMax';
            Ensemble.DivStr = 'Divergence from ensemble decision';
            OptInline1 = ['Perf' OPa 'PerfCrit && Entropy >= EntropyCrit'];
            OptInline2 = ['OrigPerf' OPb 'OptPerf || (OrigPerf == OptPerf && OrigEntropy > OptEntropy)'];
            Ensemble.OptInlineFunc1  = inline(OptInline1, 'Perf', 'PerfCrit', 'Entropy', 'EntropyCrit');
            Ensemble.OptInlineFunc2  = inline(OptInline2, 'OrigPerf', 'OptPerf', 'OrigEntropy', 'OptEntropy');

    end
end

param.SubSpaceStrategy  = SubSpaceStrategy;
param.SubSpaceCrit      = SubSpaceCrit;
param.EnsembleStrategy  = Ensemble;

if act > 0, param = nk_EnsembleStrategy2_config(param, SVM, MODEFL, [], parentstr); end