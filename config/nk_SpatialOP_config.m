function nk_SpatialOP_config(PREPROC, parentstr, stepind )

ds = 'spatial operation';
actstr=[]; actmnu=[]; 

% Get description of parameters
d = nk_GetParamDescription2(NM, PREPROC, 'SpatialOP');

if isfield(PREPROC,'SPATIAL') && ~isempty(PREPROC.SPATIAL)
    
    lact = length(PREPROC.SPATIAL);
    
    %% Display current preprocessing sequence
    fprintf('\n\n')
    cprintf('*black','SPATIAL OPERATIONS PIPELINE \n')
    cprintf('*black','=========================== ')
    for i=1:numel(d.PREPROC.spatialfiltering)
        stepstr = sprintf('Step %g: %s',i,d.PREPROC.spatialfiltering{i});
        if strcmp(PREPROC.SPATIAL{i}.cmd,'extfeat')
            cprintf([1 0.5 0],'\n%s',stepstr)
        else
            if i==stepind, 
                fprintf('\n'); cprintf('*black','=> %s',stepstr), 
            else
                fprintf('\n%s', stepstr);
            end
        end
    end
    
    cmdstr = sprintf('Add %s|Remove %s', ds, ds);
    cmdmnu = [2 3];
    [actstr, actmnu] = ConcatMenu(actstr, actmnu, cmdstr, cmdmnu);
    
    if lact > 1 
        
        if stepind == lact
            cmdstr = sprintf('Insert %s|Replace current %s|Modify current %s|<< Previous step', ds, ds, ds);
            cmdmnu = [4 5 6 7];
            [actstr, actmnu] = ConcatMenu(actstr, actmnu, cmdstr, cmdmnu);
             
        elseif stepind == 1
            cmdstr = sprintf('Insert %s|Replace current %s|Modify current %s|Next step >>',ds, ds, ds);
            cmdmnu = [4 5 6 8];
            [actstr, actmnu] = ConcatMenu(actstr, actmnu, cmdstr, cmdmnu);
        else
            cmdstr = sprintf('Insert %s|Replace %s|Modify current %s|Next step >>|<< Previous step', ds, ds, ds);
            cmdmnu = [4 5 6 8 7];
            [actstr, actmnu] = ConcatMenu(actstr, actmnu, cmdstr, cmdmnu);
        end
        cmdstr = sprintf('Change order of %s',ds);
        cmdmnu = 9;
        [actstr, actmnu] = ConcatMenu(actstr, actmnu, cmdstr, cmdmnu);
        if lact > 3
            cmdstr = 'Go to step ...'; cmdmnu = 10;
            [actstr, actmnu] = ConcatMenu(actstr, actmnu, cmdstr, cmdmnu);
        end
    else
        cmdstr = sprintf('Insert %s|Replace %s|Modify %s',ds, ds, ds);
        cmdmnu = [4 5 6];
        [actstr, actmnu] = ConcatMenu(actstr, actmnu, cmdstr, cmdmnu);
    end
    
    titlestr = sprintf('Spatial OP pipeline generator: %g / %g steps.', stepind, lact);

else
    cmdstr = sprintf('Add %s',ds);
    cmdmnu = 2;
    [actstr, actmnu] = ConcatMenu(actstr, actmnu, cmdstr, cmdmnu);
    titlestr = sprintf('Preprocessing pipeline generator: no steps defined.');
end

mestr = titlestr; navistr = [parentstr ' >>> ' mestr]; cprintf('*blue','\n\nYou are here: %s >>> ',parentstr); 
act = nk_input(titlestr, 0,'mq', actstr, actmnu);

switch act
     
    case 2 % Add Preprocessing step
        [PREPROC, EF] = config_AddReplaceModifyStep(NM, varind, PREPROC, [], 0, EF, navistr);
        if isfield(PREPROC,'ACTPARAM')
            if strcmp(PREPROC.ACTPARAM{end}.cmd,'extfeat')
                stepind = numel(PREPROC.ACTPARAM)-1;
            else
                stepind = numel(PREPROC.ACTPARAM);
            end
        end
    case 3 % Remove preprocessing step
        if lact > 1
            if strcmp(PREPROC.ACTPARAM{stepind}.cmd,'rankfeat')
                % Delete also extfeat
                remind = [stepind stepind+1];
            else
                remind = stepind;
            end
            PREPROC.ACTPARAM(remind) = [];
        else
            PREPROC = rmfield(PREPROC,'ACTPARAM');
        end
        if stepind > 1, 
            if strcmp(PREPROC.ACTPARAM{stepind-1}.cmd,'extfeat') && stepind - 2 > 0
                stepind = stepind - 2;
            else
                stepind = stepind - 1;
            end
        end
    case 4 % Insert Preprocessing step
        [PREPROC, EF] = config_AddReplaceModifyStep(NM, varind, PREPROC, stepind, 2, EF, navistr);
    
    case 5 % Replace current step
        [PREPROC, EF] = config_AddReplaceModifyStep(NM, varind, PREPROC, stepind, 1, EF, navistr);
        
    case 6 % Modify current step
        [PREPROC, EF] = config_AddReplaceModifyStep(NM, varind, PREPROC, stepind, 0, EF, navistr);
    
    case 7 % Go to previous step
        if strcmp(PREPROC.ACTPARAM{stepind-1}.cmd,'extfeat') && stepind - 2 > 0
            stepind = stepind - 2;
        else
            stepind = stepind - 1;
        end
    
    case 8 % Go to next step
        if strcmp(PREPROC.ACTPARAM{stepind}.cmd,'rankfeat') 
            if stepind + 2 <= numel(PREPROC.ACTPARAM)  
                stepind = stepind + 2;
            end
        else
            stepind = stepind + 1;
        end
        
    case 9 % Change order of preprocessing steps
        neworder = nk_input('Define new order of preprocessing steps',0,'e',[],[1, lact]);
        fl = true; 
        for i = 1:numel(neworder)
            if strcmp(PREPROC.ACTPARAM{neworder(i)}.cmd,'rankfeat') 
                if i == numel(neworder) || (i < numel(neworder) && ~strcmp(PREPROC.ACTPARAM{neworder(i+1)}.cmd,'extfeat')), fl = false;  end
            end
        end
        if fl, 
            PREPROC.ACTPARAM = PREPROC.ACTPARAM(neworder);
        else
            errordlg('Reordering cannot be performed! Weighting-based feature generation has to be executed right after the computation of the feature weighting / ranking!')
        end
        
    case 10
        tstepind = nk_input('Go to preprocessing step',0,'w1',stepind);
        if strcmp(PREPROC.ACTPARAM{tstepind}.cmd,'extfeat') 
            stepind = tstepind -1;
        else
            stepind = tstepind;
        end

end

if ~isfield(PREPROC,'BINMOD'), PREPROC.BINMOD = 1; end

% -------------------------------------------------------------------------
function [PREPROC, EF] = config_AddReplaceModifyStep(NM, varind, PREPROC, stepind, ...
                                                replflag, EF, navistr)

if ~isempty(PREPROC) && isfield(PREPROC,'ACTPARAM'), ...
    lact = length(PREPROC.ACTPARAM); 
else
    lact = 0; 
end
if ~exist('stepind','var') || isempty(stepind)
    modflag = 0; stepind = lact + 1;
else
    if replflag == 2, modflag = 0; else modflag = 1; end
end

if ~exist('replflag','var') || isempty(replflag) replflag = 0; end

actstr = []; actmnu = []; CURACT2 = [];

if modflag && ~replflag
    
    CURACT = PREPROC.ACTPARAM{stepind};
    cmdstr = CURACT.cmd;
    switch cmdstr
        case 'impute'
            cmd = 2;
        case 'correctnuis'
            cmd = 3;
        case 'scale'
            cmd = 4;
        case 'normalize'
            cmd = 5;
        case 'standardize'
            cmd = 6;
        case {'discretize','symbolize','caim_discretize','cacc_discretize'}
            cmd = 7;
        case 'reducedim'
            cmd = 8;
        case 'labelimpute'
            cmd = 9;
        case 'elimzero'
            cmd = 10;
        case 'rankfeat'
            cmd = 11;
            CURACT2 = PREPROC.ACTPARAM{stepind+1};
        case 'remmeandiff'
            cmd = 12;
        case 'unitnormalize'
            cmd = 13;
        case 'remvarcomp'
            cmd = 14;
    end
   
else    
    fld = fieldnames(EF);
    for z = 1:numel(fld)
        if EF.(fld{z})
            switch fld{z}
                 case 'impute'
                    nans = isnan(NM.Y{varind}); 
                    if sum(nans(:))
                        snan_rows = sum(nans,2); snan_mat = sum(snan_rows);
                        cmdstr = sprintf('Impute missing values (%g NaNs in %1.1f%% of cases) ',snan_mat, sum(snan_rows>0)*100/size(nans,1));    cmdmnu = 2;
                    else
                        continue
                    end
                case 'correctnuis'
                    if isfield(NM,'covars') && ~isempty(NM.covars)
                        cmdstr = 'Correct data for nuisance covariates using partial correlations'; cmdmnu = 3;
                    else
                        continue;
                    end
                case 'scale'
                    cmdstr = 'Scale data';                                                          cmdmnu = 4;
                case 'normalize'
                    if isfield(NM,'covars') && ~isempty(NM.covars)
                        cmdstr = 'Normalize to group mean';                                         cmdmnu = 5;
                    else
                        continue;
                    end
                case 'standardize'
                    cmdstr = 'Standardize data';                                                    cmdmnu = 6;
                case 'binning'
                    cmdstr = 'Apply binning method to data';                                        cmdmnu = 7;
                case 'reducedim'
                    cmdstr = 'Apply dimensionality reduction method to data';                       cmdmnu = 8;
                case 'labelimpute'
                    slnan = sum(isnan(NM.label));
                    if slnan
                        cmdstr = 'Propagate labels to unlabeled training data';                     cmdmnu = 9;
                    else
                        continue
                    end
                case 'elimzero'
                    cmdstr = 'Prune non-informative columns from data matrix';                      cmdmnu = 10;
                case 'rankfeat'
                    cmdstr = 'Rank / Weight features';                                              cmdmnu = 11;
                case 'remmeandiff'
                    if isfield(NM,'covars') && ~isempty(NM.covars)
                        cmdstr = 'Correct group offsets from global group mean';                    cmdmnu = 12;
                    else
                        continue;
                    end
                case 'unitnormalize'
                    cmdstr = 'Normalize data to unit vector';                                       cmdmnu = 13;
                case 'remvarcomp'
                    cmdstr = 'Extract variance components from data';                               cmdmnu = 14;
            end
            [actstr, actmnu] = ConcatMenu(actstr, actmnu, cmdstr, cmdmnu);
        end
    end    
    cmd = nk_input(['Select preprocessing type [Step ' ...
        num2str(stepind) ' / ' num2str(lact) ']'],0,'mq',actstr,actmnu);
    if ~cmd, return; end 
    CURACT = [];
   
end

switch cmd
    case 2
        CURACT = config_impute(NM, varind, CURACT, navistr);
    case 3
        CURACT = config_covars(NM, varind, CURACT, navistr);
    case 4
        CURACT = config_scaling(CURACT, navistr);
    case 5
        CURACT = config_groupnorm(NM, CURACT);
    case 6
        CURACT = config_standard(CURACT, navistr);
    case 7
        CURACT = config_binning(CURACT);
    case 8
        CURACT = config_dimred(CURACT, navistr);
    case 9
        CURACT = config_labelimpute(CURACT, navistr);
    case 10
        CURACT = config_elimzero(CURACT, navistr);
    case 11
        CURACT = config_rankfeat( NM, varind, CURACT, navistr );
        CURACT2 = config_waction( NM, varind, CURACT2, navistr );
    case 12
        CURACT = config_remmeandiff(NM, CURACT);
    case 13
        CURACT = config_unitnorm( CURACT, navistr );
    case 14
        CURACT = config_remvarcomp( NM, varind, CURACT, navistr );
end

switch replflag
    case 2 % Insert
        tACTPARAM = PREPROC.ACTPARAM(stepind:end);
        stepper = 1; if ~isempty(CURACT2), stepper = 2; end
        PREPROC.ACTPARAM(stepind+stepper:end+stepper) = tACTPARAM;
end
if ~isempty(CURACT),PREPROC.ACTPARAM{stepind} = CURACT; end
if ~isempty(CURACT2),PREPROC.ACTPARAM{stepind+1} = CURACT2; end

% -------------------------------------------------------------------------
function [actstr, actmnu] = ConcatMenu(actstr, actmnu, cmdstr, cmdmnu)

if isempty(actstr),
    actstr = cmdstr;
    actmnu = cmdmnu;
else
    actstr = [ actstr '|' cmdstr ];
    actmnu = [ actmnu cmdmnu ];
end
