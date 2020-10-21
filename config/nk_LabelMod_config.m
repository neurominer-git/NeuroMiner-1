function [ LABELMOD, act ] = nk_LabelMod_config( LABELMOD, parentstr, defaultsfl )

LABELMODDEF_SCALE = 1; LABELMODDEF_POLY = 1;

if ~exist('defaultsfl','var') || isempty(defaultsfl); defaultsfl = false; end

if ~defaultsfl
    
    if isfield(LABELMOD,'TARGETSCALE'), LABELMODDEF_SCALE = LABELMOD.TARGETSCALE;     end
    if isfield(LABELMOD,'POLYNOM'),     LABELMODDEF_POLY  = LABELMOD.POLYNOM;     end
    
    switch LABELMODDEF_SCALE
        case 0
            LABELMODSTR_SCALE = 'disabled'; 
        case 1
            LABELMODSTR_SCALE = 'enabled'; 
    end
    
    switch LABELMODDEF_POLY
        case 1
            LABELMODSTR_POLY = 'no transformation';
        otherwise
            LABELMODSTR_POLY = sprintf('label^%g',LABELMODDEF_POLY); 
    end
    
    menustr = ['Scale labels [ ' LABELMODSTR_SCALE ' ]|', ...
               'Exponentially transform labels [ ' LABELMODSTR_POLY ' ]'];
    menuact = 1:2;
    
    nk_PrintLogo
    mestr = 'Process-level label modification'; 
    navistr = [parentstr ' >>> ' mestr]; cprintf('*blue','\nYou are here: %s >>> ',parentstr); 
    act = nk_input(mestr,0,'mq', menustr, menuact);
    
    switch act
        case 1
            LABELMODDEF_SCALE = ~LABELMODDEF_SCALE;
        case 2
            if nk_input('Compute exponential transform of labels',0,'m','exponential transform|no transform',[1,0])
                LABELMOD.POLYNOM = nk_input('Specify exponent',0,'r',LABELMODDEF_POLY); 
            elseif isfield(LABELMOD,'POLYNOM')
                LABELMOD = rmfield(LABELMOD,'POLYNOM');
            end
    end
    
end

LABELMOD.TARGETSCALE = LABELMODDEF_SCALE;