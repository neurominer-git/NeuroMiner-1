function [CURACT, act] = nk_Unitnorm_config(CURACT, parentstr, defaultsfl)

CURACT.UNITNORM = 1; tMETHOD = 1; ZeroOut = 1; CALIBUSE = 2;

if ~exist('defaultsfl','var') || isempty(defaultsfl); defaultsfl = false; end

if ~defaultsfl
    
    if isfield(CURACT,'METHOD'), tMETHOD = CURACT.METHOD; end
    if isfield(CURACT,'zerooutflag'), ZeroOut = CURACT.zerooutflag; end
    if isfield(CURACT,'CALIBUSE'), CALIBUSE = CURACT.CALIBUSE; end
    
    switch tMETHOD
        case 1
            METHODSTR = 'L1 norm';
        case 2
            METHODSTR = 'L2 norm';
    end
    
    if ZeroOut == 1
        ZEROOUTSTR = 'yes'; 
    else
        ZEROOUTSTR = 'no';
    end
    
    menustr = ['Select normalization method [ ' METHODSTR ' ]' ...
               '|Zero-out completely non-finite features [ ' ZEROOUTSTR ' ]'];
    menuact = [1,2];
    
    [menustr, menuact] = nk_CheckCalibAvailMenu_config(menustr, menuact, CALIBUSE);
    
    nk_PrintLogo
    mestr = 'Unit-normalization'; navistr = [parentstr ' >>> ' mestr]; cprintf('*blue','\nYou are here: %s >>> ',parentstr); 
    act = nk_input(mestr,0,'mq', menustr, menuact);
    
    switch act
        case 1
            if CURACT.UNITNORM
                tMETHOD = nk_input('Select normalization method',0,'m','L1 norm|L2 norm',[1,2],tMETHOD);
            end
        case 2
             if ZeroOut == 1, ZeroOut = 2; elseif ZeroOut == 2, ZeroOut = 1; end
        case 1000
            CALIBUSE = nk_AskCalibUse_config(mestr, CALIBUSE); 
    end    
else
     act = 0;
end

CURACT.METHOD = tMETHOD;
CURACT.zerooutflag = ZeroOut;
CURACT.CLAIBUSE = CALIBUSE;