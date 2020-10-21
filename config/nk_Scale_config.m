function [ SCALE, act ] = nk_Scale_config(SCALE, parentstr, defaultsfl)

AcMatFl = 0; ZeroOne = 1; ZeroOut = 1; SCALE.xscale = 1; CALIBUSE = 2;

if ~exist('defaultsfl','var') || isempty(defaultsfl); defaultsfl = false; end

if ~defaultsfl
   
    if isfield(SCALE,'AcMatFl'), AcMatFl = SCALE.AcMatFl; end
    if isfield(SCALE,'ZeroOne'), ZeroOne = SCALE.ZeroOne; end
    if isfield(SCALE,'zerooutflag'), ZeroOut = SCALE.zerooutflag; end
    if isfield(SCALE,'CALIBUSE'), CALIBUSE = SCALE.CALIBUSE; end
    if AcMatFl
        ACMATSTR = 'Across matrix'; ACMATDEF = 1;
    else
        ACMATSTR = 'Each feature independently'; ACMATDEF = 2;
    end
    
    if ZeroOne == 1
        ZEROONESTR = '0, 1'; 
    else
        ZEROONESTR = '-1, 1';
    end
    
    if ZeroOut == 1
        ZEROOUTSTR = 'yes'; 
    else
        ZEROOUTSTR = 'no';
    end
        
    menustr = ['Scale across features or each feature independently [ ' ACMATSTR ' ]|', ...
               'Define scale range [ ' ZEROONESTR ' ]|', ...
               'Zero-out completely non-finite features [ ' ZEROOUTSTR ' ]'];
    menuact = 1:3;
    
    [menustr, menuact] = nk_CheckCalibAvailMenu_config(menustr, menuact, CALIBUSE);
    
    nk_PrintLogo
    mestr = 'Scaling'; navistr = [parentstr ' >>> ' mestr]; cprintf('*blue','\nYou are here: %s >>> ',parentstr); 
    act = nk_input(mestr,0,'mq', menustr, menuact);
    
    switch act
        case 1    
           AcMatFl = nk_input('Scale across matrix or each feature independently',0, 'm', ...
                'Across Matrix|Each Feature Independently',[1 0],ACMATDEF);
        case 2         
           ZeroOne = nk_input('Scaling range',0, 'm', 'from 0 to 1|from -1 to 1',[1 2],ZeroOne);
        case 3
           if ZeroOut == 1, ZeroOut = 2; elseif ZeroOut == 2, ZeroOut = 1; end
        case 1000
           CALIBUSE = nk_AskCalibUse_config(mestr, CALIBUSE); 
    end
else
    act = 0;
end
SCALE.AcMatFl       = AcMatFl;
SCALE.ZeroOne       = ZeroOne;
SCALE.zerooutflag   = ZeroOut;
SCALE.CALIBUSE      = CALIBUSE;