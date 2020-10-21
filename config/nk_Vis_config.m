function [VIS, act] = nk_Vis_config(VIS, PREPROC, M, defaultsfl, parentstr)
global NM 

if ~exist('defaultsfl','var') || isempty(defaultsfl), defaultsfl=0; end
if ~exist('M','var') || isempty(M), M=1; end
thresh = 0; fwhm = []; flgProj = false; 
PERM = struct('flag',0,'nperms',1000,'sigflag',0,'mode',1);
normfl = false; if any(strcmp({'LIBSVM','LIBLIN'}, NM.TrainParam.SVM.prog)), normfl = true; end
normdef = 1; if ~isfield(VIS,'norm'), VIS.norm = normdef; end   
if ~defaultsfl
    
    if ~exist('VIS','var') || isempty(VIS) , VIS = nk_Vis_config([], PREPROC, M, 1); end
    D = NM.datadescriptor{M}; 
    sourcestr = D.source; 
    menustr = []; menuact = [];
    
    if ~strcmp(sourcestr,'matrix')
        %% Gaussian smoothing setup
        if ~isfield(PREPROC,'SPATIAL') || isempty(PREPROC.SPATIAL)
            if isfield(VIS,'fwhm') && ~isempty(VIS.fwhm)
                fwhmstr = sprintf('yes (FHWM = %g)', VIS.fwhm);
                fwhmdef = VIS.fwhm; fwhmflagdef = 1;
            else
                fwhmstr = 'no';
                fwhmdef = 8; fwhmflagdef = 2;
            end
            menustr = sprintf('%sApply gaussian smoothing to weight vectors [ %s ]|', menustr, fwhmstr);
            menuact = [menuact 3];
        end
        
        %% Thresholding setup
        meanthreshstr = 'undefined'; sethreshstr = 'undefined'; threshstr = 'undefined'; threshdef = 0;
        if isfield(VIS,'thresh')
            if  VIS.thresh
                threshdef = 1; threshstr = 'yes'; threshtype = {'Percentile','Absolute Value','None'};
                if isfield(VIS,'mean'), meanthreshstr = sprintf('Type: %s; Value: %g; Logical operation: %s',threshtype{VIS.mean.type}, VIS.mean.val, VIS.mean.logop); end
                if isfield(VIS,'se'), sethreshstr = sprintf('Type: %s; Value: %g; Logical operation: %s',threshtype{VIS.se.type}, VIS.se.val, VIS.se.logop); end
            else
                threshstr = 'no';
            end
        end
        
        menustr = sprintf('%sThreshold weight vectors [ %s ]|', menustr, threshstr); menuact = [ menuact 4 ];
        if isfield(VIS,'thresh') && VIS.thresh
            menustr = sprintf('%sApply thresholding to mean weight vector images [ %s ]|', menustr, meanthreshstr); menuact = [ menuact 5 ];
            menustr = sprintf('%sApply thresholding to standard error images [ %s ]|', menustr, sethreshstr); menuact = [ menuact 6 ];
        end
    end
    
    %% Weight vector normalization setup
    if normfl
        if isfield(VIS,'norm'), normdef = VIS.norm; end
        normstr = {'yes','no'};
        menustr = sprintf('%sNormalize weight vectors [ %s ]|', menustr, normstr{normdef}); menuact = [ menuact 8 ];
    end
    
    %% Permuation setup
    permstr = 'Permutation mode disabled'; sigstr = ''; permmodestr='';
    if isfield(VIS,'PERM')
        if VIS.PERM.flag
            permstr = sprintf('yes (%g permutations)', VIS.PERM.nperms);
            if isfield(VIS.PERM,'mode')
                switch VIS.PERM.mode
                    case 1
                        permmodestr = ': labels';
                    case 2
                        permmodestr = ': features';
                    case 3
                        permmodestr = ': labels && features';
                end
            else
                VIS.PERM.mode = 3;  permmodestr = ': labels && features';
            end
       
            if iscell(PREPROC)
                zPREPROC = PREPROC{M};
            else
                zPREPROC = PREPROC;
            end
            for i=1:numel(zPREPROC.ACTPARAM)
                if strcmp(zPREPROC.ACTPARAM{i}.cmd,'reducedim'), flgProj=true; end
            end
            if flgProj
                if isfield(VIS.PERM,'sigflag') && VIS.PERM.sigflag
                    sigstr = ', back-project significant features only';
                else
                    sigstr = ', back-project all features';
                end
            end
        else
            permstr = 'no';
        end
    end
    
    menustr = sprintf('%sDerive Z scores and P values using permuatation analysis [ %s%s%s ]', menustr, permstr, permmodestr, sigstr); menuact = [ menuact 7 ];

    %% CONFIGURATION
    mestr = 'Visualization parameters'; navistr = [parentstr ' >>> ' mestr]; cprintf('*blue','\n\nYou are here: %s >>> ',parentstr); 
    nk_PrintLogo
    act = nk_input(mestr,0,'mq', menustr, menuact);

    switch act
      
        case 3
            fwhmflag = nk_input('Do you want to apply Gaussian smoothing to the weight vectors ?',0,'yes|no',[1,2], fwhmflagdef);
            if fwhmflag == 1
                VIS.fwhm = nk_input('Specify FWHM width [in mm]',0,'e', fwhmdef);
            else
                VIS.fwhm = [];
            end
        case 4
            switch NM.modeflag
                case 'classification'
                    strmode = 'discriminative';
                case 'regression'
                    strmode = 'predictive';
            end
            VIS.thresh = nk_input(['Apply thresholding to ' strmode ' images?'],0,'yes|no',[1,0], threshdef);
        case 5
            typ = []; val = []; logop = [];
            if isfield(VIS,'mean'); typ = VIS.mean.type; val = VIS.mean.val; logop = VIS.mean.logop; end
            VIS.mean = config_threshold('mean image', typ, val ,logop);
        case 6
            typ = []; val = []; logop = [];
            if isfield(VIS,'se'); typ = VIS.se.type; val = VIS.se.val; logop = VIS.se.logop; end
            VIS.se = config_threshold('standard error image', typ, val, logop);
        case 7
            VIS.PERM.flag = nk_input('Enable permutation mode',0,'yes|no',[1,0],1);
            if VIS.PERM.flag 
                VIS.PERM.nperms = nk_input('# of permutations',0,'i',100);
                if flgProj
                    VIS.PERM.sigflag = nk_input('Limit pattern generation to combinations of significant weights',0,'yes|no',[1,0],1); 
                else
                    VIS.PERM.sigflag = 0; 
                end
                VIS.PERM.mode = nk_input('Permutation mode',0,'m','Labels|Features (within-label)|Labels & Features',1:3,VIS.PERM.mode);
            end
        case 8
            if VIS.norm == 1, VIS.norm = 2; else, VIS.norm = 1; end
            
    end
else
    VIS.flipfeats = [];
    VIS.fwhm = fwhm;
    VIS.thresh = thresh;
    VIS.PERM = PERM;
    VIS.norm = normdef;
    act = 0;
end


function thresh = config_threshold(strimg, threshtype, threshval, threshlogop)

def = [];
if exist('threshtype','var') && ~isempty(threshtype)
    def = threshtype;
end

thresh.type = nk_input(['Threshold type for ' strimg],0,'m', ...
                            ['Percentile|' ...
                            'Absolute Value|' ...
                            'None'],1:3,def);

def = [];
if exist('thresval','var') && ~isempty(threshval)
    def = threshval;
end

switch thresh.type
    case 1
        thresh.val = nk_input(['Define percentile(s) for threshold of ' strimg],0,'e',def);
    case 2
        thresh.val = nk_input(['Define absolute value(s) for threshold of ' strimg],0,'e',def);
    case 3
        thresh.val = [];
end


if thresh.type == 1 || thresh.type == 2
    def = [];
    if exist('threslogop','var') && ~isempty(threshlogop) && ~threshlogop
        def = threshlogop;
    end
    switch length(thresh.val)
        case 1
            thresh.logop = nk_input('Define logical operation for thresholding', 0, 'm', ...
                ['>|' ...
                '>=|' ...
                '<|' ...
                '<='],1:4, def);
        case 2
            thresh.logop = nk_input('Define logical operation for thresholding', 0, 'm', ...
                ['>|' ...
                '>=|' ...
                '<|' ...
                '<=|' ...
                '<>|' ...
                '<=>|' ...
                '><|' ...
                '>=<'],1:8, def);
    end
else
    thresh.logop = 0;
end
