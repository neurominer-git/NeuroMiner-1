function [CURACT, act] = nk_GroupNorm_config(NM, varind, CURACT, parentstr, defaultsfl)

if ~exist('defaultsfl','var') || isempty(defaultsfl); defaultsfl = false; end
IND = []; ZeroOut = 1;
if ~defaultsfl
    
    if isfield(CURACT,'IND') && ~isempty(CURACT.IND), IND = CURACT.IND; end
    if isfield(CURACT,'zerooutflag'), ZeroOut = CURACT.zerooutflag; end
   
    if isempty(IND)
        INDSTR = 'undefined';
    else
        INDSTR = NM.covnames{IND};
    end
    
    if ZeroOut == 1, ZEROOUTSTR = 'yes'; else ZEROOUTSTR = 'no'; end
    
    menustr = ['Select group index vector for normalization (each group is labeled by a unique value in the vector) [ ' INDSTR ' ]|' ...
                'Zero-out completely non-finite features [ ' ZEROOUTSTR ' ]'];

    menuact = 1:2;
    
    nk_PrintLogo
    mestr = 'Group normalization setup'; navistr = [parentstr ' >>> ' mestr]; cprintf('*blue','\nYou are here: %s >>> ',parentstr); 
    act = nk_input(mestr,0,'mq', menustr, menuact);
    
    switch act
        case 1
            IND = nk_SelectCovariateIndex(NM,varind,1);
        case 2
            ZeroOut = nk_input('Zero-out completely non-infinite features',0, 'yes|no',[1 2],ZeroOut);
    end
else
    act = 0;
end

CURACT.IND = IND;
CURACT.zerooutflag = ZeroOut;


    
