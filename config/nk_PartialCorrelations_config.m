function [CURACT, act ] = nk_PartialCorrelations_config(NM, CURACT, varind, parentstr, defaultsfl)

if ~exist('defaultsfl','var') || isempty(defaultsfl); defaultsfl = false; end
% Defaults
METHOD_DEF      = 1;
COVAR_DEF       = 1;
MCOVARUSE_DEF   = 1;
MCOVAR_DEF      = [];
MCOVARLABEL_DEF  = 1;
COVDIR_DEF      = 1;
INTERCEPT_DEF   = 2;
BETAEXT_DEF     = [];
SUBGROUP_DEF    = [];

if ~defaultsfl
    
    % Get information from CURACT if available
    if isfield(CURACT,'METHOD'),      METHOD_DEF      = CURACT.METHOD; end
    if isfield(CURACT,'COVAR'),       COVAR_DEF       = CURACT.COVAR; end
    if isfield(CURACT,'SUBGROUP'),    SUBGROUP_DEF    = CURACT.SUBGROUP; end
    if isfield(CURACT,'INTERCEPT'),   INTERCEPT_DEF   = CURACT.INTERCEPT; end
    if isfield(CURACT,'COVDIR'),      COVDIR_DEF      = CURACT.COVDIR; end
    if isfield(CURACT,'BETAEXT'),     BETAEXT_DEF     = CURACT.BETAEXT; end
    if isfield(CURACT,'MCOVARUSE'),   MCOVARUSE_DEF   = CURACT.MCOVARUSE; end
    if isfield(CURACT,'MCOVAR'),      MCOVAR_DEF      = CURACT.MCOVAR; end
    if isfield(CURACT,'MCOVARLABEL'), MCOVARLABEL_DEF = CURACT.MCOVARLABEL; end;
    COVAR_STR = strjoin(NM.covnames(COVAR_DEF),', ');
    SUBGROUP_MNU1 = []; SUBGROUP_MNU2 = [];
    menuact = [ 1 2 ];
    
    if METHOD_DEF == 1
        % Partial Correlations
        menuact = [ menuact 3 4 ];
        METHOD_STR = 'Partial Correlations';
        if INTERCEPT_DEF == 2,          INTERCEPT_STR = 'yes'; else     INTERCEPT_STR = 'no'; end
        if COVDIR_DEF == 1,             COVDIR_STR = 'attenuate'; else  COVDIR_STR = 'increase'; end; 
    elseif METHOD_DEF==2
        % ComBat
        menuact = [ menuact 5 ]; MCOVARLABEL_MNU = []; MCOVAR_MNU = [];
        METHOD_STR = 'ComBat';
        if MCOVARUSE_DEF == 1
            MCOVARUSE_STR = 'yes';
            if isempty( MCOVAR_DEF )
                MCOVAR_STR = 'none';
            else
                MCOVAR_STR = strjoin(NM.covnames(MCOVAR_DEF),', ');
            end
            MCOVAR_MNU = sprintf('|Define retainment covariates [ %s ]', MCOVAR_STR);
            if MCOVARLABEL_DEF == 1,    MCOVARLABEL_STR = 'yes'; else, MCOVARLABEL_STR = 'no'; end
            MCOVARLABEL_MNU = sprintf('|Include NM label in variance retainment [ %s ]', MCOVARLABEL_STR);
            menuact = [ menuact 6 7 ];
        else
            MCOVARUSE_STR = 'no'; 
        end
        BETAEXT_DEF = [];
    end
    
    % Do we have to deal with external betas for partial correlations?
    if ~isempty(BETAEXT_DEF),   
        BETAEXT_STR = 'yes';
        if isfinite(BETAEXT_DEF), 
            BETAEXT_MAT = sprintf('%g x %g matrix defined', size(BETAEXT_DEF,1), size(BETAEXT_DEF,2)); 
        else
            BETAEXT_MAT = 'undefined'; 
        end
        BETAEXT_MNU = sprintf('|Define beta coeficients [ %s ]',BETAEXT_MAT);
        menuact = [ menuact 8 9 ];
        
    else
        BETAEXT_STR = 'no'; 
        BETAEXT_MNU = [];
        if METHOD_DEF == 1, menuact = [menuact 9]; end
        if ~isempty(SUBGROUP_DEF), 
            SUBGROUP_STR = 'yes'; 
            if isfinite(SUBGROUP_DEF), 
                SUBGROUP_MAT = sprintf('vector with %g case(s) defined', sum(SUBGROUP_DEF)); 
            else
                SUBGROUP_MAT = 'undefined';
            end
            SUBGROUP_MNU2 = sprintf('|Provide index to training cases for computing betas [ %s ]', SUBGROUP_MAT );
            menuact = [menuact 10 11];
        else
            SUBGROUP_STR = 'no'; SUBGROUP_MNU2 = [];
            menuact = [menuact 10];
        end
        SUBGROUP_MNU1 = sprintf('|Define subgroup of training cases for computing betas [ %s ]',  SUBGROUP_STR);
    end
    
    switch METHOD_DEF 
        case 1
            menustr = ['Select method [ ' METHOD_STR ' ]', ...
               '|Select covariates from NM covariate matrix [ ' COVAR_STR ' ]', ...
               '|Use intercept in partial correlation analysis [ ' INTERCEPT_STR ' ]', ...
               '|Attenuate or increase covariate effects [ ' COVDIR_STR ' ]' ...
               '|Use externally-computed beta coeficients [ ' BETAEXT_STR ' ]' ...
               BETAEXT_MNU ...
               SUBGROUP_MNU1 ...
               SUBGROUP_MNU2];
           
        case 2
              menustr = ['Select method [ ' METHOD_STR ' ]', ...
               '|Define vector in NM covariate matrix indexing batch effects [ ' COVAR_STR ' ]', ...
               '|Retain variance effects during correction[ ' MCOVARUSE_STR ' ]', ...
               MCOVAR_MNU ...
               MCOVARLABEL_MNU ...
               SUBGROUP_MNU1 ...
               SUBGROUP_MNU2];        
    end
           
    nk_PrintLogo
    mestr = 'Residualization setup'; navistr = [parentstr ' >>> ' mestr]; cprintf('*blue','\nYou are here: %s >>> ',parentstr); 
    act = nk_input(mestr,0,'mq', menustr, menuact);
    
    switch act
        
        case 1
            if METHOD_DEF == 1, METHOD_DEF= 2; else, METHOD_DEF= 1; end
        
        case 2
            COVAR_DEF = nk_SelectCovariateIndex(NM, COVAR_DEF, 1);
   
        case 3
            if INTERCEPT_DEF == 2, INTERCEPT_DEF = 1; elseif INTERCEPT_DEF == 1, INTERCEPT_DEF = 2; end
        
        case 4
            if COVDIR_DEF == 1, COVDIR_DEF = 2; elseif COVDIR_DEF == 2, COVDIR_DEF = 1; end
            
        case 5
             if MCOVARUSE_DEF == 1, MCOVARUSE_DEF = 2; else, MCOVARUSE_DEF = 1; end
            
        case 6
            MCOVAR_DEF = nk_SelectCovariateIndex(NM, MCOVAR_DEF, 1);
            if ~MCOVAR_DEF , MCOVAR_DEF = []; end
        
        case 7
           if MCOVARLABEL_DEF == 1, MCOVARLABEL_DEF = 2; else, MCOVARLABEL_DEF = 1; end
            
        case 8
            if ~isfield(CURACT,'BETAEXT'), 
                CURACT.BETAEXT = NaN; 
            elseif isfield(CURACT,'BETAEXT'), 
                CURACT = rmfield(CURACT,'BETAEXT'); 
            end
            
        case 9
            if INTERCEPT_DEF
                defc = numel(COVAR_DEF) + 1;
            else
                defc = numel(COVAR_DEF) ;
            end
            CURACT.BETAEXT = nk_input('Define precompute beta matrix',0,'e',[],[defc, size(NM.Y{varind},2)]);

        case 10
              if ~isfield(CURACT,'SUBGROUP'), 
                CURACT.SUBGROUP = NaN; 
            elseif isfield(CURACT,'SUBGROUP'), 
                CURACT = rmfield(CURACT,'SUBGROUP'); 
              end
              
        case 11
            CURACT.SUBGROUP = logical(nk_input('Define logical index vector to select cases for beta computation',0,'e',[],[numel(NM.label),1]));
    end
end
CURACT.METHOD       = METHOD_DEF;
CURACT.COVAR        = COVAR_DEF;
CURACT.MCOVARUSE    = MCOVARUSE_DEF;
CURACT.MCOVARLABEL  = MCOVARLABEL_DEF;
CURACT.MCOVAR       = MCOVAR_DEF;
CURACT.INTERCEPT    = INTERCEPT_DEF;
CURACT.COVDIR       = COVDIR_DEF;
