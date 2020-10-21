function [ PRUNE, PX, act ] = nk_Prune_config(PRUNE, PX, parentstr, defaultsfl)

PruneZero = 1; PruneNan = 1; PruneInf = 1; PrunePerc = [];

if ~exist('defaultsfl','var') || isempty(defaultsfl); defaultsfl = false; end

if ~defaultsfl
   
    if isfield(PRUNE,'zero'), PruneZero = PRUNE.zero; end
    if isfield(PRUNE,'nan'),  PruneNan = PRUNE.nan; end
    if isfield(PRUNE,'inf'),  PruneInf = PRUNE.inf; end
    if isfield(PRUNE,'perc') && ~isempty(PRUNE.perc), PrunePerc = PRUNE.perc; end
    
    if PruneZero == 1, PRUNEZEROSTR = 'yes'; else PRUNEZEROSTR = 'no'; end
    if PruneNan  == 1, PRUNENANSTR  = 'yes'; else PRUNENANSTR  = 'no'; end
    if PruneInf  == 1, PRUNEINFSTR  = 'yes'; else PRUNEINFSTR  = 'no'; end
    if ~isempty(PrunePerc)
        PercDef = 1;
        PRUNEPERCSTR = ['yes, ' nk_ConcatParamstr(PrunePerc)];
    else
        PercDef = 2;
        PRUNEPERCSTR = 'no';
    end
        
    menustr = ['Prune zero features [ ' PRUNEZEROSTR ' ]|', ...
               'Prune features with NaNs [ ' PRUNENANSTR ' ]|',...
               'Prune features with Infs [ ' PRUNEINFSTR ' ]|',...
               'Prune features with single-value percentage over cutoff [ ' PRUNEPERCSTR ' ]'];
    menuact = 1:4;
    
    nk_PrintLogo
    mestr = 'Pruning setup'; navistr = [parentstr ' >>> ' mestr]; cprintf('*blue','\nYou are here: %s >>> ',parentstr); 
    act = nk_input(mestr,0,'mq', menustr, menuact);
    
    switch act
        case 1
            if PruneZero == 1, PruneZero = 2; elseif PruneZero == 2, PruneZero = 1; end
        case 2         
            if PruneNan == 1, PruneNan = 2; elseif PruneNan == 2, PruneNan = 1; end
        case 3
            if PruneInf == 1, PruneInf = 2; elseif PruneInf == 2, PruneInf = 1; end
        case 4
            if PercDef == 1, PercDef = 2; elseif PercDef == 2, PercDef = 1; end
            if PercDef == 1 
                PrunePerc = nk_input('Define percentage(s) to determine single-value cutoff',0, 'i', PrunePerc);
                PX = nk_AddParam(PrunePerc, 'PrunePerc', 1, []);
            else
                PrunePerc = [];
            end
    end
else
    act = 0;
end
PRUNE.zero = PruneZero;
PRUNE.nan  = PruneNan;
PRUNE.inf  = PruneInf;
PRUNE.perc = PrunePerc;
% Generate parameter array for preprocessing pipeline runner
if exist('PX','var') && ~isempty(PX) && ...
        isfield (PX,'Px') && ...
        isfield(PX.Px,'Params') && ...
        ~isempty(PX.Px.Params)
    PX.opt = allcomb(PX.Px.Params,'matlab'); 
else
    PX.opt = [];
end
