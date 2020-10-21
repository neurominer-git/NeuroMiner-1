function [act, param] = nk_Filter_config(act, param, SVM, MODEFL, MULTI, defaultsfl, parentstr)
global NM
%
% Defaults:
% ---------
Filter.flag             = 0; % Activate feature filtering
Filter.MinNum           = 1; % Minimum # of features
% Default ranking method
switch MODEFL 
    case 'classification'
        Filter.type     = 3;
        Filter.Pearson  = 1;
    case 'regression'
        Filter.type     = 2;
end
Filter.binmode          = 1;
Filter                  = nk_CostFun_config(Filter, SVM, MODEFL, 1); % Default cost function
% Default subspace strategy
Filter.SubSpaceFlag     = 1;
Filter.RankThresh       = 95;
Filter                  = nk_SubSpaceStrategy_config(Filter, SVM, MODEFL, 1);
Filter                  = nk_EnsembleStrategy2_config(Filter, SVM, MODEFL, 1);
limfeat                 = 2;
Filter.SubSpaceStepping = 0;
Filter.PFE              = nk_ProbalisticFea_config([],1);

if ~exist('defaultsfl','var') || isempty(defaultsfl),  defaultsfl = 0; end;

% -------------------------------------------------------------------------
if ~isempty(act) || ~defaultsfl
    
    % -------------------------------------------------------------------------
    % Retrieve current values from param
    if isfield(param,'Filter'), 
        Filter = param.Filter; 
    end
         
    % Get Parameter descriptions
    d = nk_GetParamDescription2([],param,'FeatFlt');
    
    % Defaults at the prompt
    if ~Filter.flag,                 flagdef = 0; else              flagdef = Filter.flag;                  end
    if ~Filter.SubSpaceFlag,         subspaceflagdef = 2; else      subspaceflagdef = Filter.SubSpaceFlag;  end
    if ~Filter.binmode,              binmodedef = 2; else           binmodedef = Filter.binmode;            end
    if ~Filter.SubSpaceStepping,     steppingdef = 2; else          steppingdef = 1;                        end
         
    menustr = sprintf('Train filter methods on CV1 partitions [ %s ]', d.FeatFltFlag); menuact = 1;
    
    if strcmp(MODEFL,'classification')
        menustr = sprintf('%s|Use algorithm output scores or label predictions [ %s ]', menustr, d.FeatMetric);                 menuact = [ menuact 2 ];
    else
        Filter.EnsembleStrategy.Metric = 1;
    end
    
    if Filter.flag 
        %menustr = sprintf('%s|Define filter evaluation mode [ %s ]', menustr, d.FeatFltSubSpaces);                              menuact = [ menuact 3 ];
        if Filter.SubSpaceFlag
            if strcmp(MODEFL,'classification') && ~isempty(MULTI) && MULTI.flag
                menustr = sprintf('%s|Optimize for multi-class settings [ %s ]', menustr, d.FeatFltBinmode);                    menuact = [ menuact 4 ];
            end
             menustr = sprintf('%s|Specify filter type [ %s ]', menustr, d.FeatFlt);                                            menuact = [ menuact 5 ];
             menustr = sprintf('%s|Define minimum # of features to be selected [ %1.0f ]', menustr, Filter.MinNum);             menuact = [ menuact 6 ];
             menustr = sprintf('%s|Specify subspace optimization strategy [ %s ]', menustr, d.FeatFltEnsStrat);                 menuact = [ menuact 7 ];
             menustr = sprintf('%s|Define target population for computing optimization parameter [ %s ]', menustr, d.FeatPop);  menuact = [ menuact 8 ];
             menustr = sprintf('%s|Define subspace stepping [ %s ]', menustr, d.FeatFltStepping);                               menuact = [ menuact 9 ];
        else
             menustr = sprintf('%s|Define percentile cut-off threshold for ranked features [ %s ]', menustr, d.FeatFltThresh);  menuact = [ menuact 10 ];
             menustr = sprintf('%s|Define filter type [ %s ]', menustr, d.FeatFlt);                                             menuact = [ menuact 5 ];
             menustr = sprintf('%s|Define minimum # of features to be selected [ %1.0f ]', menustr, Filter.MinNum);             menuact = [ menuact 6 ];
        end
    else
        if isfield(param,'Wrapper') && isfield(param.Wrapper,'SubSpaceFlag') && param.Wrapper.SubSpaceFlag
            param.Wrapper.SubSpaceFlag = 0;
            param.Wrapper.EnsembleStrategy.type = 0;
            cprintf('red','Ensemble-based subspace optimization in Wrapper module DISABLED');
        end
    end
    nk_PrintLogo
    mestr = 'Filter-based model selection setup'; navistr = sprintf('%s\n\t>>> %s',parentstr, mestr); cprintf('*blue','\nYou are here: %s >>> ',parentstr);
    
    act = nk_input(mestr, 0,'mq', menustr, menuact);
    
    switch act
        
        case 1
            Filter.flag = ~flagdef;
            
        case 2
            Filter.EnsembleStrategy.Metric = nk_input('Use predicted labels or algorithm scores for ensemble construction',0,'m', ...
                                                        ['Predicted labels (Hard decision ensemble)|' ...
                                                         'Algorithm scores (Soft decision ensemble)'], 1:2, Filter.EnsembleStrategy.Metric);
        case 3
            Filter.SubSpaceFlag = nk_input('Filter evaluation mode', 0,'Subspace|Ranking',[1,0], subspaceflagdef);
            if Filter.SubSpaceFlag, 
                Filter = nk_SubSpaceStrategy_config(Filter, SVM, MODEFL, [], navistr); 
            end
        case 4
            Filter.binmode = nk_input('Feature ranking mode',0, 'binary|multi-group',[1,0], binmodedef);
            
        case 5
            switch MODEFL
                case 'classification'
                     Filter.type = nk_input('Filter method',0,'m', ...
                        ['Mutual information-based (FEAST toolbox)|' ...
                        'Maximum Relevance Minimum Redundancy (Peng et al)|' ...
                        'Pearson/Spearman correlation coefficient|' ...
                        'Simba (margin based method)|' ...
                        'Gflip (margin based method)|' ...
                        'IMRelief (margin-based method)|' ...
                        'Increasing feature subspace cardinality|' ...
                        'Random subspace sampling|'...
                        'Relief (margin-based method for classification)|' ...
                        'FScore|' ...
                        'Bhattacharya distance'],[1:5,7,9,11,12,13,14], ...
                        Filter.type);
                    
                case 'regression'
                     Filter.type = nk_input('Filter method',0,'m', ...
                        ['Mutual information-based (FEAST toolbox)|' ...
                        'Maximum Relevance Minimum Redundancy (Peng et al)|' ...
                        'Pearson/Spearman correlation coefficient|' ...
                        'RGS|' ...
                        'Increasing feature subspace cardinality|' ...
                        'Random subspace sampling|'...
                        'Relief (margin-based method for regression)'],[1 2 3 10 9 11 12], ...
                        Filter.type); 
            end
            switch Filter.type
                case 1
                    acti=1; while acti > 0, [acti, Filter] = nk_FEAST_config( Filter, [], navistr); end
                case 2
                    acti=1; while acti > 0, [acti, Filter] = nk_MRMR_config( Filter, [], navistr ); end
                case 3
                    if isfield(Filter,'Pearson'), fltdef = Filter.Pearson; else, fltdef = 2; end
                    Filter.Pearson = nk_input('Select type univariate correlation analysis',0,'Pearson|Spearman',[1,2],fltdef);
                case 4
                    Filter = nk_Simba_config(Filter, 0, 0, length(unique(NM.label)));
                case 5
                    Filter = nk_Gflip_config(Filter, 0 ,0, length(unique(NM.label)));
                    if Filter.SubSpaceFlag
                        promptstr = 'Limit features subspaces to features selected by G-flip?';
                    else
                        promptstr = 'Use only features selected by G-flip?';
                    end
                    Filter.gflip.limfeat = nk_input(promptstr,0, 'yes|no',[1,0],limfeat);
                case 7
                    Filter = nk_SetupIMReliefParams_config(NM, Filter);
                case 10
                    Filter = nk_RGS_config(Filter);
                case 11
                    Filter = nk_RSS_config(Filter);
                case 12
                    Filter.Relief.k = nk_input('Define number of nearest neigbours for RELIEF',0,'e',10);
            end
        case 6
            Filter.MinNum = nk_input('Minimum # of features to use',0,'e',Filter.MinNum);
        case 7
            Filter = nk_EnsembleStrategy2_config(Filter, SVM, MODEFL, [], navistr);
        case 8
            Filter = nk_CostFun_config(Filter, SVM, MODEFL);
        case 9
            SubSpaceSteppingFlag = nk_input('Define feature subspace stepping',0,'yes|no',[1,0], steppingdef);
            if SubSpaceSteppingFlag
                Filter.SubSpaceStepping = nk_input('Percentage of maximum feature number per subspace',0,'w1',Filter.SubSpaceStepping);
            else
                Filter.SubSpaceStepping = 0;
            end
        case 10
            Filter.RankThresh = nk_input('Define percentile cut-off threshold for ranked features',0,'e',Filter.RankThresh);
    end
end
param.Filter = Filter;

