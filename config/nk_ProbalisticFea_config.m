function [PFE, act] = nk_ProbalisticFea_config(PFE, parentstr, defaultsfl)

if ~exist('defaultsfl','var') || isempty(defaultsfl),  defaultsfl = 0; end;

tPFE.flag            = 0;   % Perform PFE
tPFE.Mode            = 1;    % Perform PFE
tPFE.Perc            = 75;   % Percentage of CV1 feature mask agreement
tPFE.TolWin          = 25;   % Percentage of CV1 feature mask agreement
tPFE.MinNum          = 1;    % Minimum # features required across CV1
tPFE.PruneFlag        = 1;    % Flag for different operations using the Cross-C1 agreement vector

if ~defaultsfl && ~isempty(PFE)
    
    if isfield(PFE,'flag'), tPFE.flag       = PFE.flag; end
    if isfield(PFE,'Mode'), tPFE.Mode       = PFE.Mode; end
    if isfield(PFE,'Perc'), tPFE.Perc       = PFE.Perc; end
    if isfield(PFE,'TolWin'), tPFE.TolWin   = PFE.TolWin; end
    if isfield(PFE,'MinNum'), tPFE.MinNum   = PFE.MinNum; end
    if isfield(PFE,'PruneFlag'), tPFE.PruneFlag = PFE.PruneFlag; end
    if isfield(PFE,'flag'), tPFE.flag = PFE.flag; end
    if ~tPFE.PruneFlag, tPFE.PruneFlag = 2; else, tPFE.PruneFlag = 1; end
    
    if tPFE.flag==1, PROBFEASTR_FLAG = 'yes'; else PROBFEASTR_FLAG = 'no'; end
    
    switch tPFE.Mode
        case 1
            PROBFEASTR_MODE = '% Cross-CV1 feature selection agreement';
        case 2
            PROBFEASTR_MODE = sprintf('%g of most consistently selected features', 100-tPFE.Perc);
        case 3
            PROBFEASTR_MODE = sprintf('%g%% of most consistently selected features',100-tPFE.Perc);
    end
    
    switch tPFE.PruneFlag
        case 1
            PROBFEASTR_PRUNEFLAG = 'Prune unselected features from each model and retrain models';
        case 2
            PROBFEASTR_PRUNEFLAG = 'Retrain models using single optimized feature mask';
    end
    
    menustr = ['Optimize feature selection across CV1 partitions [ ' PROBFEASTR_FLAG ' ]']; menuact = 1;
    
    if tPFE.flag 
        menustr = [menustr '|Probabilistic feature extraction mode [ ' PROBFEASTR_MODE ' ]']; menuact = [menuact 2];
        menustr = [menustr '|Appy consistency-based ranking to [ ' PROBFEASTR_PRUNEFLAG ' ]']; menuact = [menuact 3];
    end
    
    nk_PrintLogo
    mestr = 'Probabilistic feature extraction setup'; navistr = [parentstr ' >>> ' mestr]; cprintf('*blue','\nYou are here: %s >>> ',parentstr); 
    act = nk_input(mestr,0,'mq', menustr, menuact);
    
    switch act
    
        case 1
            tPFE.flag = ~tPFE.flag;
        case 2
            if tPFE.flag == 1
                Mode  = nk_input('Cross-CV1 feature selection mode',0,'mq', ...
                                            ['Cross-CV1 feature selection agreement with tolerance window|' ...
                                            'Absolute number of most consistently selected features|' ...
                                            'Percentage cutoff for most consistently selected features'],1:3,tPFE.Mode);
                if Mode >0, tPFE.Mode=Mode; else return; end
                switch tPFE.Mode
                    case 1
                        tPFE.Perc    = nk_input('Define minimum percentage of cross-CV1 agreement as cutoff for feature selection',0,'i', tPFE.Perc);
                        tPFE.TolWin  = nk_input('Define percentage tolerance window to be tested if this threshold does not return any features',0,'i',tPFE.TolWin);
                        tPFE.MinNum  = nk_input('Define minimum no. of features',0,'i',tPFE.MinNum);
                    case 2
                        tPFE.Perc    = nk_input('Define number of features to be selected from consistency-based ranking list',0,'i', tPFE.Perc);
                    case 3
                        tPFE.Perc    = 100 - nk_input('Define cutoff [%] for consistency-based feature selection ',0,'i', 100-tPFE.Perc);
                end

            end    
        case 3
            tPFE.PruneFlag = nk_input('Define how the cross-CV1 feature agreement data should be used',0, 'm', ...
                ['Retrain all CV1 models after pruning them from unselected features|',...
                 'Retrain all CV1 models using the same selected feature space'], [1,0], tPFE.PruneFlag);
    end
end

PFE = tPFE;
