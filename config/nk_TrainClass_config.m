% ============================================================================
% FORMAT res = nk_TrainClass_config(act, varind, parentstr)
% ============================================================================
% Interface for defining the parameters of the training and prediction process
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Nikolaos Koutsouleris, 05/2018

function [act, varind] = nk_TrainClass_config(act, varind, parentstr)
global NM

if ~exist('act','var'), act = []; end
menustr = []; menuact = []; 

%% Select variate to work on
if ~exist('varind','var') || isempty(varind), 
    varind = 1; 
    if isfield(NM,'TrainParam') && ...
            isfield(NM.TrainParam,'FUSION') && isfield(NM.TrainParam.FUSION,'M')
        varind = NM.TrainParam.FUSION.M(1);
    end
end

%% Check whether TrainParams already exist
if ~isfield(NM,'TrainParam') 
    % Create default NM parameters space
    nk_CVpartition_config(true);
    NM.TrainParam.STACKING.flag = 2; 
    NM.TrainParam.FUSION.flag   = 0;
    NM.TrainParam.FUSION.M      = 1;
    NM.TrainParam.SVM           = nk_LIBSVM_config(NM,[],1);
    NM.TrainParam.SVM.prog      = 'LIBSVM';
    NM.TrainParam.SVM           = nk_Kernel_config(NM.TrainParam.SVM,1);
    NM.TrainParam.SVM.GridParam = 1;
    if strcmp(NM.modeflag, 'regression'), NM.TrainParam.SVM.GridParam = 18; end
    NM.TrainParam.MULTI.flag    = 0;
    NM.TrainParam               = nk_Grid_config(NM.TrainParam, NM.TrainParam.SVM, true);
    [~,NM.TrainParam.RFE]       = nk_RFE_config([], NM.TrainParam, NM.TrainParam.SVM, NM.modeflag, NM.TrainParam.MULTI, NM.TrainParam.GRD, 1); 
    NM.TrainParam.verbosity     = 1;
end

% Check whether the number of modalities specified in FUSION.M matches 
% the number of modalities in NM (Users were copying their TrainParam 
% structures from one NM workspace to the other and in cases where a data
% fusion was implemented in on NM workspace but was not possible in another
% workspace NM crashed.
if isfield(NM,'TrainParam') && ...
            isfield(NM.TrainParam,'FUSION') && isfield(NM.TrainParam.FUSION,'M')
    nM1 = numel(NM.Y);
    nM2 = numel(NM.TrainParam.FUSION.M);
    if nM1 < nM2 || nM1 < max(NM.TrainParam.FUSION.M) 
        NM.TrainParam.FUSION.M = 1;
        NM.TrainParam.FUSION.flag = 0;
    end
end

% for compatibility reasons
if ~isfield(NM.TrainParam,'STACKING') || ~NM.TrainParam.STACKING.flag
    NM.TrainParam.STACKING.flag=2; 
end

% Adjust varind when in stacking mode and remove any SPATIAL operations
% from current modality
if NM.TrainParam.STACKING.flag == 1
    NM.TrainParam.FUSION.M = 1;
    NM.TrainParam.FUSION.flag=0;
    varind = 1;
    if isfield(NM.TrainParam.PREPROC{varind},'SPATIAL')
        NM.TrainParam.PREPROC{varind} = rmfield(NM.TrainParam.PREPROC{varind},'SPATIAL');
    end
end

%% Create further default configurations
if ~isfield(NM.TrainParam,'PREPROC'),
    % Create PREPROC structure
    nan_in_pred=false;          if sum(isnan(NM.Y{varind}(:)))>0, nan_in_pred=true; end
    nan_in_label=false;         if sum(isnan(NM.label(:)))>0, nan_in_label=true; end
    NM.TrainParam.PREPROC{1}    = DefPREPROC(NM.modeflag,nan_in_pred,nan_in_label);
    NM.TrainParam.VIS{1}        = nk_Vis_config([], NM.TrainParam.PREPROC, 1, 1); 
else
    switch NM.TrainParam.FUSION.flag
        case 3 % If late fusion has been activated create STRAT structure
            nY = numel(NM.Y);
            for i=1:nY
                if isempty(NM.TrainParam.STRAT{i})
                    NM.TrainParam.STRAT{i}.SVM          = NM.TrainParam.SVM;
                    NM.TrainParam.STRAT{i}.GRD          = NM.TrainParam.GRD;
                    if numel(NM.TrainParam.PREPROC) == i
                        NM.TrainParam.STRAT{i}.PREPROC  = NM.TrainParam.PREPROC{i};
                        NM.TrainParam.STRAT{i}.VIS      = NM.TrainParam.VIS{i};
                    else
                        NM.TrainParam.STRAT{i}.PREPROC  = NM.TrainParam.PREPROC{1};
                        NM.TrainParam.STRAT{i}.VIS      = NM.TrainParam.VIS{1};
                    end
                    NM.TrainParam.STRAT{i}.RFE          = NM.TrainParam.RFE;
                    NM.TrainParam.STRAT{i}.MULTI        = NM.TrainParam.MULTI;
                end
            end
        otherwise
            nY = numel(NM.Y);
            nP = numel(NM.TrainParam.PREPROC); 
            for i=1:nY;
                if i>nP
                    NM.TrainParam.PREPROC{i} = DefPREPROC(NM.modeflag);
                    NM.TrainParam.VIS{i}    = nk_Vis_config([], NM.TrainParam.PREPROC{i}, i, 1); 
                end
            end
    end
end

%% Check data entry status
if isfield(NM.TrainParam,'FUSION') && NM.TrainParam.FUSION.flag == 3
    STATUS = nk_CheckFieldStatus(NM,{'TrainParam','cv'},{'RAND', 'SAV', 'OOCV', 'META', 'STACKING'});
    STATUS = nk_CheckFieldStatus(NM.TrainParam.STRAT{varind},{'PREPROC','SVM','GRD','RFE','MULTI','VIS'}, [], [], STATUS);
else
    STATUS = nk_CheckFieldStatus(NM,{'TrainParam','cv'},{'STACKING','RAND','PREPROC','SVM','GRD','RFE','MULTI','VIS','SAV','OOCV'});
end
switch STATUS.PREPROC
    case '...'
        STATUS.FEATGEN = '...';
    otherwise 
        STATUS.FEATGEN = '???';
end

switch NM.TrainParam.verbosity
    case 0
        verbostr = 'No output';
    case 1
        verbostr = 'Detailed output';
end

if ~exist('act','var') || isempty(act)
    
    %% Check whether there are analyses that have been completed and make stacking options available
    s = nk_GetNMStatus(NM);
    if ~isempty(s.completed_analyses) && any(s.completed_analyses & s.isequalcv_analyses) && sum(s.nmodal_analyses)>1
        menustr = [ menustr sprintf('Define meta-learning/stacking options [ %s ]|', STATUS.STACKING) ]; menuact = [ menuact 18 ]; 
    end
    
    %% Check whether more than one variate are available and make data fusion options available
    if length(NM.Y)>1 && NM.TrainParam.STACKING.flag == 2
        % Make data fusion option available
        if isfield(NM.TrainParam,'FUSION')
            fusemode = NM.TrainParam.FUSION.flag;
            switch fusemode
                case 0 
                    fusstr = ['[ Disabled. All operations will be applied to modality #' num2str(varind) ' ]'];
                case 1
                    fusstr  = '[ BEFORE preprocessing => ';
                case 2
                    fusstr  = '[ AFTER preprocessing => ';
                case 3
                    fusstr  = '[ Decision-based fusion (bagging) => ';
            end
            if fusemode         
                fusesel = NM.TrainParam.FUSION.M;
                nF = numel(fusesel);
                if nF > 3
                    fusstr = [ fusstr sprintf('%g modalities ]',nF) ];
                else
                    fusstr = [ fusstr 'modalities:' ];
                    for i=1:nF, fusstr = sprintf('%s #%g +', fusstr, NM.TrainParam.FUSION.M(i)); end
                    fusstr = [fusstr(1:end-2) ' ]'];
                end
            end
        else
            fusstr = '[ undefined ]';
            fusemode = 0;
        end
        fusstr = ['Define data fusion options ' fusstr '|']; menustr = [menustr fusstr]; menuact = [menuact 1];
        switch fusemode 
            case {2,3}
                % Make active modality selection option available
                descstr = [ ' (' NM.datadescriptor{varind}.desc ')']; 
                varstr = ['Set active modality for configuration [ #' num2str(varind) descstr ' ]|']; menuact = [menuact 2];
                if fusemode == 3, NM.TrainParam.META.flag = 1; end
                resetstr = '';
                if fusemode == 3,
                    resetstr = 'Reset all modality configurations to current modality setup|';  menuact = [ menuact 4];
                end
                menustr = [menustr varstr resetstr];
        end
    else
        fusemode=0;
    end
     
    % Check for regression / classification experiment and number of groups
    multistr = ''; multiflag = false;
    if strcmp(NM.modeflag,'classification') 
        classtr = ['Classification algorithm [ ' STATUS.SVM ' ]|'];
        if numel(unique(NM.label(~isnan(NM.label)))) > 2
            multistr = ['Multi-class settings [ ' STATUS.MULTI ' ]|'];
            multiflag = true;
        else
            NM.TrainParam.MULTI.flag=0;
        end
    else
        classtr = ['Prediction algorithm [ ' STATUS.SVM ' ]|'];
        multistr = ''; NM.TrainParam.MULTI.flag=0;
    end
   
    if isfield(NM,'OOCV') && ~isempty(NM.OOCV)
        oocvstr = ['OOCV prediction parameters [ ' STATUS.OOCV ' ]|'];
        oocvflag = true;
    else
        oocvstr='';
        oocvflag = false;
    end

    menustr = [menustr 'Cross-validation settings [ ' STATUS.cv ' ]|']; menuact = [menuact 3];
    
    SVM = [];
    if isfield(NM.TrainParam,'FUSION') && NM.TrainParam.FUSION.flag == 3
        flFUS = true;
        flSVM = isfield(NM.TrainParam.STRAT{varind},'SVM');
        flGRD = isfield(NM.TrainParam.STRAT{varind},'GRD');
        flPREPROC = isfield(NM.TrainParam.STRAT{varind},'PREPROC');
        if flSVM, SVM = NM.TrainParam.STRAT{varind}.SVM; end
    else
        flFUS = false;
        flSVM = isfield(NM.TrainParam,'SVM');
        flGRD = isfield(NM.TrainParam,'GRD');
        flPREPROC = isfield(NM.TrainParam,'PREPROC');
        if flSVM, SVM = NM.TrainParam.SVM; end
    end
    
    flx = flSVM && flGRD && flPREPROC;
    
    menustr = [ menustr 'Preprocessing pipeline [ ' STATUS.PREPROC ' ]|' classtr ]; menuact = [ menuact 5:6 ];
    
    if flx
        
        if (any(strcmp(SVM.prog,{'MikRVM','MKLRVM','MVTRVR','MVTRVM','GLMFIT'})) && any(strcmp(SVM.kernel.kernstr,{' -t 0','lin','linear'}))) || ...
                ( strcmp(SVM.prog,'matLRN') && ( ~isfield(NM.TrainParam.GRD.matLearn,'Params') || isempty(NM.TrainParam.GRD.matLearn.Params)) ), 
            if flFUS
                NM.TrainParam.STRAT{varind}.GRD.PX = []; NM.TrainParam.STRAT{varind}.GRD.n_params = 0;
            else
                NM.TrainParam.GRD.PX = []; NM.TrainParam.GRD.n_params = 0;
            end
        else
            menustr = [ menustr 'Learning algorithm parameters [ ' STATUS.GRD ' ]|' ]; menuact = [ menuact 7 ]; 
        end
        % For the predictive sequence optimizer we don't need any ensemble
        % learning optimization.
        if ~strcmp(SVM.prog,'SEQOPT')
            menustr = [ menustr ...
                    'Ensemble generation strategies [ ' STATUS.RFE ' ]|'];
            menuact = [menuact 8];
        end
        if multiflag 
            menustr = [ menustr multistr];
            menuact = [ menuact 9 ]; 
        end
    end   
    menustr = [ menustr ... 
                'Visualization options [ ' STATUS.VIS ' ]|' ...
                'Model saving options [ ' STATUS.SAV ' ]|' ...
                oocvstr ...
                'Define verbosity level [ ' verbostr ' ]|' ...
                'Inspect workspace|' ...
                'Save parameter template|' ...
                'Load training template'];
                
                
    menuact = [ menuact 11 12 ];
    if oocvflag, menuact = [ menuact 13 ]; end
    menuact = [ menuact 16 19 14 15 ];
    
    nk_PrintLogo
    mestr = 'Define parameter template'; navistr = sprintf('%s\n\t>>> %s',parentstr, mestr); cprintf('*blue','\nYou are here: %s >>> ',parentstr);
    if fusemode > 1, fprintf('\n'); cprintf('*blue','==> CONFIGURATION OF MODALITY #%g: %s', varind, descstr ); end
    act = nk_input(mestr, 0, 'mq', menustr, menuact);
end

switch act
    
    case 0
        return
    case 1
        if isfield(NM.TrainParam,'FUSION'), 
            fusedef = NM.TrainParam.FUSION.flag ;
        else
            fusedef = 1; 
        end
        NM.TrainParam.FUSION.flag = nk_input('Select multimodal fusion strategy',0,'m', ...
            ['No fusion|' ...
             'Early fusion -> Modality concatenation BEFORE feature preprocessing|' ...
             'Intermediate fusion -> Modality concatenation AFTER preprocessing|' ...
             'Late fusion -> Decision-based data fusion (bagging)'],0:3, fusedef); 
        if NM.TrainParam.FUSION.flag 
            NM.TrainParam.FUSION.M = nk_SelectVariateIndex(NM, 0, 1, 1 );
            if numel(NM.TrainParam.FUSION.M) == 1
                NM.TrainParam.FUSION.flag = false;                
            end
            if isempty(find(NM.TrainParam.FUSION.M == varind)),
                varind = NM.TrainParam.FUSION.M(1);
            end
            if NM.TrainParam.FUSION.flag == 3
                if ~isfield(NM.TrainParam,'STRAT')
                    NM.TrainParam.STRAT = cell(numel(NM.Y),1);
                end
            end
        else
            % Only one modality has to be selected for analysis
            NM.TrainParam.FUSION.M = nk_SelectVariateIndex(NM, 1, 1, 1 );
            varind = NM.TrainParam.FUSION.M;
        end
        
    case 2
        if isfield(NM.TrainParam,'FUSION') && NM.TrainParam.FUSION.flag ;
            M = NM.TrainParam.FUSION.M ;
        else
            M = [];
        end
        varind = nk_SelectVariateIndex(NM, 1, 1, 1, M);
    
    case 3
        act = 1; while act>0, act = nk_CVpartition_config; end
    
    case 4
        fl = questdlg('Are you sure you want to overwrite all modality configuration with the current modality setup?',mestr,'Yes','No','No');
        if strcmp(fl,'Yes')
            STRAT = NM.TrainParam.STRAT(varind);
            NM.TrainParam.STRAT(NM.TrainParam.FUSION.M) = STRAT;
        end
            
    % PREPROCESSING ================================================================================================================================================= 
    case 5
        if ~isfield(NM,'TrainParam') || ...
                ~isfield(NM.TrainParam,'PREPROC') || ...
                varind > numel(NM.TrainParam.PREPROC)
            NM.TrainParam.PREPROC{varind}.FEATSEL.active=0;
            NM.TrainParam.PREPROC{varind}.BINMOD=1;
        end
        if isfield(NM.TrainParam,'FUSION') && NM.TrainParam.FUSION.flag == 3
            if ~isfield(NM.TrainParam.STRAT{varind},'PREPROC'), NM.TrainParam.STRAT{varind}.PREPROC = []; end
            act = 1; stepind = 1; while act>0, [NM.TrainParam.STRAT{varind}.PREPROC, act, stepind] = nk_Preproc_config(NM.TrainParam.STRAT{varind}.PREPROC, varind, navistr, stepind); end
        else
            if isfield(NM.TrainParam,'FUSION') && NM.TrainParam.FUSION.flag == 1
                varind = NM.TrainParam.FUSION.M(1); %Always use the first modality to setup the parameters for the concatenated feature space
            end
            if ~isfield(NM.TrainParam,'PREPROC'), NM.TrainParam.PREPROC{varind} = []; end
            act = 1; stepind = 1; while act>0, [NM.TrainParam.PREPROC{varind}, act, stepind] = nk_Preproc_config(NM.TrainParam.PREPROC{varind}, varind, navistr, stepind); end
        end
        
    % ML ALGORITHM SELECTION =========================================================================================================================================    
    case 6
        if isfield(NM.TrainParam,'FUSION') && NM.TrainParam.FUSION.flag == 3
            if ~isfield(NM.TrainParam.STRAT{varind},'SVM'), NM.TrainParam.STRAT{varind}.SVM = []; end
             act = 1; while act>0, ...
                     [NM.TrainParam.STRAT{varind}.SVM, act] = nk_Model_config(NM.TrainParam.STRAT{varind}.SVM, NM.TrainParam.STRAT{varind}, navistr, varind); end
        else
            if ~isfield(NM.TrainParam,'SVM'), NM.TrainParam.SVM = []; end
            act = 1; while act>0, [NM.TrainParam.SVM, act ] = nk_Model_config(NM.TrainParam.SVM, NM.TrainParam, navistr, varind); end
        end
        
    % ML OPTIMIZATION STRATEGIES =====================================================================================================================================      
    case 7
        if isfield(NM.TrainParam,'FUSION') && NM.TrainParam.FUSION.flag == 3
            if ~isfield(NM.TrainParam.STRAT{varind},'GRD'), NM.TrainParam.STRAT{varind} = nk_Grid_config(NM.TrainParam.STRAT{varind}, NM.TrainParam.STRAT{varind}.SVM, true); end
            act = 1; while act>0, ...
                    [ NM.TrainParam.STRAT{varind}, act ] = nk_Grid_config(NM.TrainParam.STRAT{varind}, NM.TrainParam.STRAT{varind}.SVM, [], navistr); end
        else
            if ~isfield(NM.TrainParam,'GRD'), NM.TrainParam = nk_Grid_config([], NM.TrainParam.SVM, true); end
            act = 1; while act>0, [ NM.TrainParam, act ] = nk_Grid_config(NM.TrainParam, NM.TrainParam.SVM, [], navistr); end
        end
    % FEATURE SELECTION ==============================================================================================x=================================================  
    case 8
        if isfield(NM.TrainParam,'FUSION') && NM.TrainParam.FUSION.flag == 3
            if ~isfield(NM.TrainParam.STRAT{varind},'RFE'), 
                [~, NM.TrainParam.STRAT{varind}.RFE ] = ...
                    nk_RFE_config([], NM.TrainParam.STRAT{varind}, ...
                        NM.TrainParam.STRAT{varind}.SVM, ...
                        NM.modeflag, ...
                        NM.TrainParam.STRAT{varind}.MULTI, ...
                        NM.TrainParam.STRAT{varind}.GRD, 1); 
            end
            acti = 1; while acti>0
                [ acti, NM.TrainParam.STRAT{varind}.RFE ] = ...
                    nk_RFE_config(acti, NM.TrainParam.STRAT{varind}, ...
                        NM.TrainParam.STRAT{varind}.SVM, ...
                        NM.modeflag, ...
                        NM.TrainParam.STRAT{varind}.MULTI, ...
                        NM.TrainParam.STRAT{varind}.GRD, [], navistr); 
            end
        else
            if ~isfield(NM.TrainParam,'RFE'), [~, NM.TrainParam.STRAT{varind}.RFE ] = nk_RFE_config([], NM.TrainParam, NM.TrainParam.SVM, NM.modeflag, NM.TrainParam.MULTI, NM.TrainParam.GRD, 1); end
            acti = 1; while acti>0, [ acti, NM.TrainParam.RFE ] = nk_RFE_config(acti, NM.TrainParam, NM.TrainParam.SVM, NM.modeflag, NM.TrainParam.MULTI, NM.TrainParam.GRD, [], navistr); end
        end
        
    % MULTI-GROUP SETTINGS =============================================================================================================================================   
    case 9
        if isfield(NM.TrainParam,'FUSION') && NM.TrainParam.FUSION.flag == 3
            if ~isfield(NM.TrainParam.STRAT{varind},'MULTI'), NM.TrainParam.STRAT{varind}.MULTI = nk_Multi_config([], true); end
            act = 1; while act>0, [ NM.TrainParam.STRAT{varind}.MULTI, act ] = nk_Multi_config(NM.TrainParam.STRAT{varind}.MULTI, [], navistr); end
        else
            if ~isfield(NM.TrainParam,'MULTI'), NM.TrainParam.MULTI = nk_Multi_config([], true); end
            act = 1; while act>0, [ NM.TrainParam.MULTI, act ] = nk_Multi_config(NM.TrainParam.MULTI,[], navistr); end
        end
        
    % DATA FUSION-BASED ENSEMBLE STRATEGIES ============================================================================================================================   
    case 10
        %nk_Ensemble_config;
        
    % VISUALIZATION ====================================================================================================================================================    
    case 11
        if isfield(NM.TrainParam,'FUSION') && NM.TrainParam.FUSION.flag == 3
            if ~isfield(NM.TrainParam,'VIS'), NM.TrainParam.STRAT{varind}.VIS = nk_Vis_config(NM.TrainParam.STRAT{varind}.VIS, NM.TrainParam.STRAT{varind}.PREPROC, 1, 1, navistr); end
            act = 1; while act>0, [ NM.TrainParam.STRAT{varind}.VIS, act] = nk_Vis_config(NM.TrainParam.STRAT{varind}.VIS, NM.TrainParam.STRAT{varind}.PREPROC, 1, [], navistr); end
        else
            if ~isfield(NM.TrainParam,'VIS'), NM.TrainParam.VIS{varind} = nk_Vis_config(NM.TrainParam.VIS{varind}, NM.TrainParam.PREPROC{varind}, varind, 1, navistr); end
            act = 1; while act>0, [ NM.TrainParam.VIS{varind}, act ] = nk_Vis_config(NM.TrainParam.VIS{varind}, NM.TrainParam.PREPROC{varind}, varind , [], navistr); end
        end 
        
    % SAVING OPTIONS ====================================================================================================================================================     
    case 12
        NM.TrainParam.SAV.savemodel = nk_input('Save models?',0,'yes|no',[1,0],2);
        NM.TrainParam.SAV.matname = nk_input('Define prefix of output MAT-files',0,'s'); 
    
    case 13
        NM.TrainParam = nk_OOCV_config(NM.TrainParam);
    
    case 14
        matname = nk_input('Filename (prefix is TRAIN)',0,'s');
        matname = ['TRAIN_' matname];
        if isfield(NM,'TrainParam')
            TrainParam = NM.TrainParam;
            save(matname,'TrainParam');
        end
        if isfield(NM,'cv'), 
            cv = NM.cv;
            save(matname,'cv','-append');
        end        
        
    case 15
        
        fl = nk_input('Overwrite current settings?',0,'yes|no',[1,0],1);
        if fl
            menuvec = []; menustr =[];
            if isfield(NM,'analysis'), 
                menustr = 'Take parameters from analysis in current NM structure|';
                menuvec = 1; 
            end
            menustr = [menustr 'Load parameters from parameter file|'];
            menuvec(end+1) = 2; 
            menustr = [menustr 'Load parameters from NM structure file'];
            menuvec(end+1) = 3; 
            act = nk_input('Source of parameters to be loaded', 0, 'm', menustr, menuvec);
                
            switch act 
                case {2,3}
                    mess = []; 
                    switch act
                        case 2
                            matname = nk_FileSelector(1,'matrix','select TRAIN mat','TRAIN.*\mat');
                            if exist(matname,'file'),load(matname); end
                        case 3
                            matname = nk_FileSelector(1,'matrix','Select NM structure file','.*\mat');
                            if exist(matname,'file'),
                                [~, matfile] = fileparts(matname);
                                fprintf('\nLoading %s as temporary structure',matfile) 
                                load(matname,'NM','TrainParam');
                                load(matname,'NM','cv');
                            end
                    end
                    if exist('TrainParam','var'), 
                        if isfield(NM,'TrainParam') && isfield(NM.TrainParam,'RAND')
                            RAND = NM.TrainParam.RAND;
                        end
                        NM.TrainParam = TrainParam; mess{1} = 'TrainParam'; 
                        if exist('RAND','var'), NM.TrainParam.RAND = RAND; end
                    end
                    if exist('cv','var') && ~isfield(NM,'cv') 
                        NM.cv = cv; mess{1+end} = 'CV'; 
                    else
                        mess{1+end} = 'CV already existed, not overwritten!'; 
                    end
                case 1
                    analind = nk_SelectAnalysis(NM); mess=[];
                    if isfield(NM.analysis{analind}.params,'TrainParam')
                        mess{1} = 'TrainParam';
                        NM.TrainParam = NM.analysis{analind}.params.TrainParam;
                    end
                    if isfield(NM.analysis{analind}.params,'cv')
                        NM.cv = NM.analysis{analind}.params.cv;
                        mess{1+end} = 'CV'; 
                    end
                    NM.cv = NM.analysis{analind}.params.cv;
            end
            if ~isempty(mess)>0, h = msgbox(mess,'Loaded parameters:','none'); end
        end
        
    case 16
        NM.TrainParam.verbosity = ~NM.TrainParam.verbosity;
    
    case 17
        act = 1; stepind = 1; while act>0, [NM.TrainParam.META, act, stepind] = nk_Ensemble_config(NM.TrainParam.STRAT{varind}.PREPROC, varind, navistr, stepind); end
    
         % META-LEARNING (STACKING) =================================================================================================================================================     
    case 18
        if ~isfield(NM.TrainParam,'STACKING'), NM.TrainParam.STACKING.flag = 2; end
        mess=[];act = 1; while act>0, [NM.TrainParam.STACKING, act, mess] = nk_Stacking_config(NM.TrainParam.STACKING, s, mess, navistr); end
        
    case 19
        nk_PrintWs(NM, NM.TrainParam)
       
end

act = 1;

function PREPROC = DefPREPROC(modeflag, nan_in_pred, nan_in_label)

if ~exist('nan_in_pred','var') || isempty(nan_in_pred), nan_in_pred = false; end
if ~exist('nan_in_label','var') || isempty(nan_in_label), nan_in_label = false; end

    PREPROC.BINMOD = 1; PREPROC.FEATSEL.active = 0;
    
    PREPROC.ACTPARAM{1}.SCALE   = nk_Scale_config([],[],1);
    PREPROC.ACTPARAM{1}.cmd     = 'scale';
    
    PREPROC.ACTPARAM{2}.PRUNE   = nk_Prune_config([],[],[],1);
    PREPROC.ACTPARAM{2}.cmd     = 'elimzero';
    
    if nan_in_pred % Adjust defaults to NaN in the predictor data
        PREPROC.ACTPARAM{1}.SCALE.zerooutflag   = 1;
        PREPROC.ACTPARAM{2}.PRUNE.nan           = 2;
        PREPROC.ACTPARAM{3}.IMPUTE              = nk_Impute_config([],[],[],[],1);
        PREPROC.ACTPARAM{3}.cmd                 = 'impute';
    end
    
    if nan_in_label % Adjust defaults to NaN in the labels data
        PREPROC.LABELMOD.LABELIMPUTE            = nk_Impute_config([], [], [], [], 1);
        PREPROC.LABELMOD.cmd                    = 'labelimpute';
    end
    
    if strcmp(modeflag,'regression'), PREPROC.LABELMOD.TARGETSCALE = true; end
    