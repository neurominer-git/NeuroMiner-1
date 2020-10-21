function [RANK, PX] = nk_Rank_config(RANK, PX, NM, varind, defaultsfl, parentstr)

if ~exist('defaultsfl','var') || isempty(defaultsfl),  defaultsfl = 0; end;
ranktype            = 1;
weightmethod        = 1;
label               = NM.label;
labeldesc           = 'NM target label';
algostr             = 'pearson';
act = 0;
vartype             = NM.datadescriptor{varind}.type;

if ~defaultsfl

    if ~isfield(RANK,'ranktype'),       RANK.ranktype = ranktype; end
    if ~isfield(RANK,'weightmethod'),   RANK.weightmethod = weightmethod; end
    if ~isfield(RANK,'label');          RANK.label = label; RANK.labeldesc = labeldesc;end
    if ~isfield(RANK,'algostr'),        RANK.algostr = algostr; RANK.Pearson = 1; end
    if ~exist('PX','var'),              PX = []; end
    if RANK.weightmethod == 1,          weightstr = 'upweight features'; else weightstr = 'downweight features'; end
    
    nk_PrintLogo
    
    mestr = 'Rank / Weight features'; navistr = [parentstr ' >>> ' mestr]; cprintf('*blue','\nYou are here: %s >>> ',parentstr); 
    if strcmp(RANK.algostr,'extern')
        act = nk_input(mestr,0,'mq', ...
            [sprintf('Choose algorithm and specify its parameters [ %s ]|', RANK.algostr) ...
             sprintf('Up- or downweight predictive features [ %s ]', weightstr) ], [1 3]);
    else
        switch RANK.algostr
            case 'pearson'
                if RANK.Pearson ==1, algostr = 'Pearson'; else, algostr = 'Spearman'; end
            otherwise
                algostr = RANK.algostr;
        end
        switch RANK.algostr
            case 'pls'
                act = nk_input(mestr,0,'mq', ...
                    [sprintf('Choose algorithm and specify its parameters [ %s ]|', algostr), ...
                     sprintf('Define PLS parameters ...')],[1 2]);
            otherwise
                act = nk_input(mestr,0,'mq', ...
                    [sprintf('Choose algorithm and specify its parameters [ %s ]|', algostr)  ...
                     sprintf('Define targets for ranking/weighting [ %s ]|', RANK.labeldesc) ...
                     sprintf('Up- or downweight predictive features [ %s ]', weightstr)],1:3);
        end
    end
    
    switch act
        case 1
            if exist('PX','var'), PX = []; end
            
            %% Ranking algorithm list depending on framework
            rstr_cl = ['Pearson correlation|' ...
                       'Spearman correlation|' ...
                       'G-flip|' ...
                       'Simba|' ... %'Large-Margin Nearest Neighbour Weighting (LMNN)|' ...
                       'IMRelief|' ...
                       'Relief for classification|' ...
                       'Linear SVC (LIBSVM)|' ...
                       'Linear SVC (LIBLINEAR)|' ...
                       'F-Score|' ...
                       'AUC|' ...
                       'FEAST|' ...
                       'ANOVA|' ...
                       'PLS|' ...
                       'External ranking'];
            rsel_cl = {'pearson','spearman','gflip','simba','imrelief','relief','libsvm','liblin','fscore','auc','feast','anova','pls','extern'};
            
            rstr_reg = ['Pearson correlation|' ...
                        'Spearman correlation|' ...
                        'RGS|' ...
                        'Linear SVR (LIBSVM)|' ...
                        'Linear SVR (LIBLINEAR)|' ...
                        'Relief for regression|' ...
                        'ANOVA|', ...
                        'PLS|', ...
                        'External ranking'];
            rsel_reg = {'pearson','spearman','rgs','libsvm','liblin','relief','anova','pls','extern'};
    
            switch RANK.ranktype
                case 1
                    switch NM.modeflag
                        case 'classification'
                            menuact = rstr_cl; menusel = rsel_cl;
                        case 'regression'
                            menuact = rstr_reg; menusel = rsel_reg;
                    end
                case 2
                    menuact = rstr_cl; menusel = rsel_cl;
                case 3
                    menuact = rstr_reg; menusel = rsel_reg;
            end

            RANK.algostr = char(nk_input('Select ranking algorithm',0,'m', menuact, menusel));
           
            switch RANK.algostr
                
                case 'gflip'
                    ngroups = numel(unique(RANK.label));
                    if isfield(RANK,'gflip') 
                        gflip = RANK.gflip; 
                    else
                        gflip = nk_Gflip_config([], 1, [], ngroups);
                    end
                    gflip = nk_Gflip_config(gflip, 0, [], ngroups); RANK.gflip = gflip;
                    if isfield(gflip.gflip,'extra_param') && ...
                        isnumeric(gflip.gflip.extra_param.beta),                 PX = nk_AddParam(gflip.gflip.extra_param.beta, 'Beta', 1, PX, 'replace');
                    end
                    
                case 'simba'
                    ngroups = numel(unique(RANK.label));
                    if isfield(RANK,'simba') 
                        simba = RANK.simba; 
                    else
                        simba = nk_Simba_config([], 1,[],ngroups);
                    end
                    simba = nk_Simba_config(simba, 0,[],ngroups); RANK.simba = simba;
                    if isfield(simba.simba,'extra_param') && ...
                            isnumeric(simba.simba.extra_param.beta),             PX = nk_AddParam(simba.simba.extra_param.beta, 'Beta',1, PX, 'replace');
                    end

                case 'idetect'
                    if isfield(RANK,'idetect') 
                        idetect = RANK.idetect; setupfl = false;
                    else
                        idetect = []; setupfl = true; 
                    end
                    RANK.idetect = nk_IMRelief_config(idetect, [], setupfl); 
                    PX = nk_AddParam(RANK.idetect.sigma, 'Sigma',1 ,PX, 'replace'); PX = nk_AddParam(RANK.idetect.lambda, 'Lambda',1, PX); 
                    
                case 'imrelief'
                    if isfield(RANK,'imrelief') 
                        imrelief = RANK.imrelief; setupfl = false;
                    else
                        imrelief = []; setupfl = true; 
                    end
                    RANK.imrelief = nk_IMRelief_config(imrelief, [], setupfl); 
                    PX = nk_AddParam(RANK.imrelief.sigma, 'Sigma',1 ,PX, 'replace'); PX = nk_AddParam(RANK.imrelief.lambda, 'Lambda',1, PX); 
                    
                case {'libsvm','liblin'}
                    if ~isfield(RANK,'SVM'), RANK.SVM = []; end;
                    RANK.SVM.prog = upper(RANK.algostr); 
                    switch RANK.ranktype
                        case 1
                            RANK.SVM.modeflag = NM.modeflag;
                            RANK.SVM.evalfunc = 8;
                        case 2
                            RANK.SVM.modeflag = 'classification'; 
                            RANK.SVM.evalfunc = 14;
                        case 3
                            RANK.SVM.modeflag = 'regression';
                            RANK.SVM.evalfunc = 18;
                    end
                    % Only linear kernel is permitted
                    RANK.SVM.kernel.kerndef = 1;
                    switch RANK.algostr
                        case 'libsvm'
                            RANK.SVM = nk_LIBSVM_config(RANK.SVM, RANK.SVM, [],[], navistr);
                            RANK.SVM.kernel.kernstr = ' -t 0';
                        case 'liblin'
                            RANK.SVM = nk_LIBLIN_config(RANK.SVM, RANK.SVM, [], navistr);
                            RANK.SVM.kernel.kernstr = 'lin';
                    end
                    ctype = nk_GetLIBSVMClassType(RANK.SVM);
                    % Check if regression parameters are needed 
                    rtype = nk_GetLIBSVMRegrType(RANK.SVM);
                    switch rtype
                        case 1
                            if isfield(RANK.SVM,'EpsParam'), EpsParam = RANK.SVM.EpsParam; else EpsParam = 0.1; end
                            RANK.SVM.EpsParam = nk_input('Define Epsilon parameter(s)',0,'e',EpsParam);     PX = nk_AddParam(RANK.SVM.EpsParam, 'Epsilon', 1, PX, 'replace');
                            
                        case 2
                            if isfield(RANK.SVM,'NuParam'), NuParam = RANK.SVM.NuParam; else NuParam = 0.5; end
                            RANK.SVM.NuParam = nk_input('Define Nu parameter(s)',0,'e',NuParam);            PX = nk_AddParam(RANK.SVM.NuParam, 'Nu', 1, PX, 'replace');
                            
                    end
                    % This is the slack / nu-SVC parameter of the SVM
                    if ~ctype || rtype, Slackdesc = 'Slack'; else Slackdesc = 'Nu (SVC)';end
                    if isfield(RANK.SVM,'SlackParam'), SlackParam = RANK.SVM.SlackParam; else SlackParam = 1; end
                    RANK.SVM.SlackParam = nk_input(['Define ' Slackdesc ' parameter(s)'],0,'e', SlackParam);PX = nk_AddParam(RANK.SVM.SlackParam, Slackdesc, 1, PX);
                   
                case 'rgs'
                    if ~isfield(RANK,'RGS'), RANK = nk_RGS_config(RANK, true); end
                    RANK = nk_RGS_config(RANK);                                                     
                    if isfield(RANK.RGS.extra_param,'k') && isnumeric(RANK.RGS.extra_param.k), 
                        PX = nk_AddParam(RANK.RGS.extra_param.k, 'K', 1, PX, 'replace'); 
                    else
                        PX = nk_AddParam([], [], [], [], 'reset');
                    end
                    if isfield(RANK.RGS.extra_param,'beta') && isnumeric(RANK.RGS.extra_param.beta), 
                        PX = nk_AddParam(RANK.RGS.extra_param.beta, 'Beta', 1, PX); 
                    end
                case 'feast'
                    act=1; while act > 0, [act, RANK ] = nk_FEAST_config(RANK, 0, navistr); end
                case 'relief'
                    RANK.Relief.k = nk_input('Define number(s) of nearest neigbours for RELIEF',0,'i',10,1);PX = nk_AddParam(RANK.Relief.k, 'K', 1, PX); 
                case 'pls'
                    if ~isfield(RANK,'PLS'), [~, RANK.PLS ] = nk_PLS_config(NM, [], [], true); end
                case 'extern' 
                    if vartype
                        if isfield(NM.datadescriptor{varind},'Yw')
                            readimg = nk_input('Do you want to read-in external ranking map',0,'mq', ...
                                ['Use weighting image in NM|' ...
                                 'Import external weighting map from file|' ...
                                 'Read-in weighting vector from MATLAB workspace'],[2,1,0],1);
                        else
                             readimg = nk_input('Do you want to read-in external ranking map',0,'mq', ...
                                ['Import external weighting map from file|' ...
                                 'Read-in weighting vector from MATLAB workspace'],[1,0],1);
                        end
                    else
                        readimg = false;
                    end
                    switch readimg
                        case 0
                            RANK.EXTERN = nk_input('Define external map/vector',0,'e',[],[1 size(NM.Y{varind},2)]);
                        case 1
                            % Currently only NIFTI/Analyze supported
                            imgtype = NM.datadescriptor{varind}.input_settings.datasource;
                            RANK.V  = [];
                            [RANK.F, RANK.V{1}] = nk_FileSelector(1,imgtype,'Select Ranking Map','.*\.nii$|.*\.img$');
                            switch imgtype
                                case {'spm','nifti'}
                                     Thresh = NM.datadescriptor{varind}.input_settings.Thresh;
                                     Vm = spm_vol(NM.brainmask{varind}); 
                                     [RANK.EXTERN, RANK.THRESH] = nk_ReturnSubSpaces(RANK.V, Vm, 1, 1, Thresh);
                                case 'surf'
                                     RANK.EXTERN = MRIread(RANK.F);
                            end
                        case 2
                            nPw = size(NM.datadescriptor{varind}.input_settings.Pw,1);
                            if nPw>1
                                mn_sel = strjoin(cellstr(NM.datadescriptor{varind}.input_settings.Pw),'|');
                                mn_act = 1:nPw;
                                RANK.V = nk_input('Select weighting image',0,'m',mn_act,mn_sel,1);
                            else
                                RANK.V = 1;
                            end
                            if isfield(RANK,'EXTERN'), RANK = rmfield(RANK,'EXTERN'); end
                            if isfield(RANK,'F'), RANK = rmfield(RANK,'F'); end
                    end
            end
            
        case 2
            switch RANK.algostr
                
                case 'pls'
                    PX = nk_AddParam([], [], [], PX,'reset');
                    if ~isfield(RANK,'PLS'), [~, RANK.PLS ] = nk_PLS_config(NM, [], [], true); end
                    act = 1; while act >0, [act, RANK.PLS, PX ] = nk_PLS_config(NM, RANK.PLS, PX, false, navistr); end  
                    RANK.label = RANK.PLS.V;
                    switch RANK.PLS.uselabel
                        case 1
                            RANK.labeldesc = 'NM target label';
                        case 2
                            RANK.labeldesc = 'NM covariate matrix';
                        case 3
                            RANK.labeldesc = 'NM user-defined matrix';
                    end
                otherwise
                    
                    %% Label definition
                    if isfield(RANK,'label')
                        flg = nk_input('Target label vector found',0,'Use existing|Define new one',[0,1],1);
                    else
                        flg = true;
                    end

                    if flg, 
                        RANK.ranktype = nk_input('Define target labels',0,'m', ...
                                            ['NM target labels|'...
                                             'Categorical data|' ...
                                             'Continuous data'],[1,2,3],ranktype);

                        if RANK.ranktype > 1
                            switch RANK.algostr
                                case 'anova'
                                    RANK.label = nk_input('Define design matrix for ANOVA',0,'e',[],[numel(NM.label),Inf]); lbstr = 'design matrix';
                                otherwise
                                    RANK.label = nk_input('Define target label vector',0,'e',[],[numel(NM.label),1]); lbstr = 'label vector';
                            end
                            RANK.labeldesc = nk_input(['Give a short description of the ' lbstr ],0,'s');
                        else
                            RANK.label = NM.label;
                            RANK.labeldesc = 'NM target label';
                        end
                    end
                    
            end
            
            if RANK.ranktype == 3
                flg = nk_input('Weight feature using only one specific subgroup?',0,'yes|no',[1,0],0);
                if flg, 
                    RANK.glabel = nk_input('Define index vector for subgroup idenitification (ones = use / zeros = not use)',0,'e',[],[numel(NM.label),1]);
                else
                    if isfield(RANK,'glabel'), RANK= rmfield(RANK,'glabel'); end
                end
            else
                if isfield(RANK,'glabel'), RANK= rmfield(RANK,'glabel'); end
            end

        case 3
            RANK.weightmethod = nk_input('Select weighting method',0, 'm', ...
                'Upweight relevant features|Downweight relevant features',[1,2], weightmethod);
    end
else
    RANK.ranktype       = ranktype;
    RANK.weightmethod   = weightmethod;
    RANK.label          = NM.label;
    RANK.labeldesc      = 'NM target label';
    RANK.rankmethod     = 8;
    RANK.algostr        = 'pearson';
end
if act, [RANK, PX] = nk_Rank_config(RANK, PX, NM, varind, [], parentstr); end
