function [TrainParam, act] = nk_Grid_config(TrainParam, SVM, defaultsfl, parentstr)
global NM

OptimFlag = 1;
if ~exist('defaultsfl','var') || isempty(defaultsfl), defaultsfl=0; end
[~, CompStr] = nk_ReturnEvalOperator(SVM.GridParam);

switch OptimFlag
    
    case 1
        
        % Define defaults
        switch NM.modeflag
            case 'classification'
                Cdefs                   = 2.^(-6:1:4);
            case 'regression'
                Cdefs                   = 2.^(-6:1:0);
        end
        Gdefs                           = 2.^(-15:2:3);
        Epsdefs                         = [ 0.05 0.075 0.1 0.125 0.15 0.2];
        Nudefs                          = [ 0.2 0.5 0.7];
        Toldefs                         = [1e-6 1e-5 1e-4 1e-3 1e-2];
        Tdefs                           = [.2 .5 .7];
        PolyCoefdefs                    = 0;
        PolyDegrdefs                    = 3;
        Neurondefs                      = [25 50 75 100];
        Leafdefs                        = logspace(1,2,10);
        Treedefs                        = [25 50 75 100 150 200];
        NumDdefs                        = 50;
        Kdefs                           = 7;
        Weightdefs                      = 1;
        CoxCutoffsdefs                  = [20:10:80];
        SEQOPTstepsdefs                 = 10;
        SEQOPTlimsLdefs                 = 50;
        SEQOPTlimsUdefs                 = 50;
        OptRegul.flag                   = 0;
        GridMaxType                     = 1;
        OptRegul.type                   = 1;
        OptRegul.lambda                 = 1;
        OptRegul.big_gamma              = .5;
        OptRegul.RegulTypeComplexity    = nk_RegulFunc_config(1);
        OptRegul.RegulTypeDiversity     = nk_RegulFunc_config(1);
        NodeSelect.mode                 = 1;
        switch CompStr
            case 'above'
                 NodeSelect.perc        = 95;
                 NodeSelect.percvec     = 75:5:95;
            case 'below'
                 NodeSelect.perc        = 5;
                 NodeSelect.percvec     = 5:5:25;
        end
        act                             = 0;
        PX                              = [];
        GRD                             = [];
        
        if ~defaultsfl 
            
            if isfield(TrainParam,'GRD'), GRD = TrainParam.GRD; end
            
            %% Get current values 
            %==============================================================
            if isfield(GRD,'Cparams'),                  Cdefs   = GRD.Cparams; end
            if isfield(GRD,'Epsparams'),                Epsdefs = GRD.Epsparams; end
            if isfield(GRD,'Nuparams'),                 Nudefs  = GRD.Nuparams; end
            if isfield(GRD,'Gparams'),                  Gdefs = GRD.Gparams; end
            if isfield(GRD,'PolyCoefparams'),           PolyCoefdefs = GRD.PolyCoefparams; end
            if isfield(GRD,'PolyDegrparams'),           PolyDegrdefs = GRD.PolyDegrparams; end
            if isfield(GRD,'Neuronparams'),             Neurondefs = GRD.Neuronparams; end
            if isfield(GRD,'Kparams'),                  Kdefs = GRD.Kparams; end
            if isfield(GRD,'Tolparams'),                Toldefs = GRD.Tolparams; end
            if isfield(GRD,'Leafparams'),               Leafdefs = GRD.Leafparams; end
            if isfield(GRD,'Treeparams'),               Treedefs = GRD.Treeparams; end
            if isfield(GRD,'NumDparams'),               NumDdefs = GRD.NumDparams; end
            if isfield(GRD,'Weightparams'),             Weightdefs = GRD.Weightparams; end
            if isfield(GRD,'CutOffparams'),             SEQOPTstepsdefs = GRD.CutOffparams; end
            if isfield(GRD,'LimsLparams'),              SEQOPTlimsLdefs = GRD.LimsLparams; end
            if isfield(GRD,'LimsUparams'),              SEQOPTlimsUdefs = GRD.LimsUparams; end
            if isfield(GRD,'CoxCutoffparams')           CoxCutoffsdefs = GRD.CoxCutoffparams; end
            if isfield(GRD,'OptRegul'),                 OptRegul = GRD.OptRegul; end
            if isfield(GRD,'NodeSelect'),               NodeSelect = GRD.NodeSelect; end
            
            %============================================================== 
            menustr = []; menuact = []; n_pars = [];
            
            %% Slack setup
            switch SVM.prog
                
                case {'LIBSVM','LIBLIN','SVMLIT','SVMPRF', 'MEXELM'}
                    ctype = nk_GetLIBSVMClassType(SVM);
                    switch ctype
                        case 1
                            Rparstr = 'Nu-SVC parameter(s)'; [Nustr, n_pars(end+1)] = nk_ConcatParamstr(Nudefs); 
                            PX = nk_AddParam(Nudefs, ['ML-' Rparstr], 2, PX);
                            menustr = sprintf('Define %s [ %s ]', Rparstr, Nustr);                                  menuact = [ menuact 3 ]; 
                        otherwise
                            Cparstr = 'Slack/Regularization parameter(s)'; [Cstr, n_pars(end+1)] = nk_ConcatParamstr(Cdefs); 
                            PX = nk_AddParam(Cdefs, ['ML-' Cparstr], 2, PX);
                            menustr = sprintf('Define %s [ %s ]', Cparstr, Cstr);                                   menuact = [ menuact 1 ]; 
                    end
                    rtype = nk_GetLIBSVMRegrType(SVM);
                    switch rtype
                        case 1
                            Rparstr = 'Eps-SVR parameter(s)'; [Epsstr, n_pars(end+1)] = nk_ConcatParamstr(Epsdefs);
                            PX = nk_AddParam(Epsdefs, ['ML-' Rparstr], 2, PX);
                            menustr = sprintf('%s|Define %s [ %s ]', menustr, Rparstr, Epsstr);                     menuact = [ menuact 2 ];
                        case 2
                            Rparstr = 'Nu-SVR parameter(s)'; [Nustr, n_pars(end+1)] = nk_ConcatParamstr(Nudefs);
                            PX = nk_AddParam(Nudefs, ['ML-' Rparstr], 2, PX);
                            menustr = sprintf('%s|Define %s [ %s ]', menustr, Rparstr, Nustr);                      menuact = [ menuact 3 ];
                    end
                    
                case 'IMRELF'
                    Cparstr = 'Lambda parameter(s)'; [Cstr, n_pars(end+1)] = nk_ConcatParamstr( Cdefs );
                    PX = nk_AddParam(Cdefs, ['ML-' Cparstr], 2, PX);
                    menustr = sprintf('Define %s [ %s ]',Cparstr, Cstr);                                            menuact = [ menuact 1 ];
                    
                case 'matLRN'
                    if ~isfield(GRD,'matLearn') || ~isfield(GRD.matLearn,'Params')
                        GRD.matLearn = nk_matLearn_config(SVM.matLRN,SVM.matLRN.learner.framework,3);
                    end
                    if ~isempty(GRD.matLearn.Params)
                        for i=1:numel(GRD.matLearn.Params)
                            parstr = GRD.matLearn.Params(i).name; [~, n_pars(end+1)] = nk_ConcatParamstr(  GRD.matLearn.Params(i).range );
                            PX = nk_AddParam(GRD.matLearn.Params(i).range, ['ML-' parstr], 2, PX);
                        end
                        str = numel(GRD.matLearn.Params);
                        menustr = sprintf('%s|Define matLearn parameters [ %g main parameters defined ]', menustr, str); 
                        menuact = [ menuact 20 ];
                    else
                        return
                    end
                case {'GLMNET','GRDBST','ROBSVM'}
                    if isfield(GRD,SVM.prog) && ~isempty(GRD.(SVM.prog).Params)
                        for i=1:numel(GRD.(SVM.prog).Params)
                            parstr = GRD.(SVM.prog).Params(i).name; [~, n_pars(end+1)] = nk_ConcatParamstr(  GRD.(SVM.prog).Params(i).range );
                            PX = nk_AddParam(GRD.(SVM.prog).Params(i).range, ['ML-' parstr], 2, PX);
                        end
                        str = sprintf('%g parameters defined', numel(GRD.(SVM.prog).Params));
                    else
                        str = 'undefined'; 
                    end
                    menustr = sprintf('%s|Define %s parameters [ %s ]|', menustr, SVM.prog, str); 
                    menuact = [ menuact 21 ];
                
                    
            end
            
            %% Kernel setup
            switch SVM.prog
                
                case 'IMRELF'
                     Gparstr = 'Sigma parameter(s)'; [Gstr, n_pars(end+1)] = nk_ConcatParamstr(Gdefs);
                     PX = nk_AddParam(Gdefs, ['ML-' Gparstr], 2, PX);
                     menustr = sprintf('%s|Define %s [ %s ]', menustr, Gparstr, Gstr);                              menuact = [ menuact 4 ];                     
                    
                otherwise
                     switch SVM.kernel.kernstr
                        case {' -t 1', ...
                            ' -t 2', ...
                            ' --t 2', ...
                            ' -t 3', ...
                            ' -t 6', ...
                            ' -t 7', ...
                            'poly', ...
                            'polynomial', ...
                            'Polynomial', ...
                            'polyN', ...
                            'hpolyN', ...
                            'rbf', ...
                            'expo', ...
                            'laplace', ...
                            'cauchy', ...
                            'cubic', ...
                            'tps', ...
                            'r', ...
                            'gauss', ...
                            'gaussian'}
                            Gparstr = 'Kernel parameter(s)'; [Gstr, n_pars(end+1)] = nk_ConcatParamstr(Gdefs);
                            PX = nk_AddParam(Gdefs, ['ML-' Gparstr], 2, PX);
                            menustr = sprintf('%s|Define %s [ %s ]', menustr, Gparstr, Gstr);                           menuact = [ menuact 4 ];
                            % Additional polynomial kernel params
                            switch SVM.kernel.kernstr
                                
                                case {' -t 1', 'poly', 'polynomial', 'Polynomial', 'hpolyN', 'polyN'}
                                    Pcparstr = 'Polynomial coefficients'; [Pcstr, n_pars(end+1)] = nk_ConcatParamstr(PolyCoefdefs);
                                    PX = nk_AddParam(PolyCoefdefs, ['ML-' Pcparstr], 2, PX);
                                    menustr = sprintf('%s|Define %s [ %s ]', menustr, Pcparstr, Pcstr);                 menuact = [ menuact 5 ];
                                    Pdparstr = 'Polynomial degree'; [Pdstr, n_pars(end+1)] = nk_ConcatParamstr(PolyDegrdefs);
                                    PX = nk_AddParam(PolyDegrdefs, ['ML-' Pdparstr], 2, PX);
                                    menustr = sprintf('%s|Define %s [ %s ]', menustr, Pdparstr, Pdstr);                 menuact = [ menuact 6 ];
                                    
                                case ' -t 3'
                                    Pcparstr = 'Sigmoid coefficients'; [Pcstr, n_pars(end+1)] = nk_ConcatParamstr(PolyCoefdefs);
                                    PX = nk_AddParam(PolyCoefdefs, ['ML-' Pcparstr], 2, PX);
                                    menustr = sprintf('%s|Define %s [ %s ]', menustr, Pcparstr, Pcstr);                 menuact = [ menuact 5 ];
                            end    
                     end             
            end
            
            %% Other param setup
            switch SVM.prog
                case {'LIBSVM','LIBLIN'}
                    if SVM.(SVM.prog).Weighting
                        Weightparstr = 'Weighting exponents'; [Weightstr, n_pars(end+1)] = nk_ConcatParamstr(Weightdefs);
                        PX = nk_AddParam(Weightdefs, ['ML-' Weightparstr], 2, PX);
                        menustr = sprintf('%s|Define %s [ %s ]', menustr, Weightparstr, Weightstr);                     menuact = [ menuact 15 ];
                    end
                case 'MEXELM'
                    Neuronparstr = 'Hidden neurons'; [Neuronstr, n_pars(end+1)] = nk_ConcatParamstr(Neurondefs);
                    PX = nk_AddParam(Neurondefs, ['ML-' Neuronparstr], 2, PX);
                    menustr = sprintf('%s|Define %s [ %s ]', menustr, Neuronparstr, Neuronstr);                         menuact = [ menuact 12 ];
                case 'kNNMEX'
                    Kparstr = 'Nearest neighbors'; [Kstr, n_pars(end+1)] = nk_ConcatParamstr(Kdefs);
                    PX = nk_AddParam(Kdefs, ['ML-' Kparstr], 2, PX);
                    menustr = sprintf('%s|Define %s [ %s ]', menustr, Kparstr, Kstr);                                   menuact = [ menuact 13 ];
                case 'BLOREG'
                    Tparstr = 'Tolerance parameter(s)'; [Tstr, n_pars(end+1)] = nk_ConcatParamstr(Toldefs);
                    PX = nk_AddParam(Toldefs, ['ML-' Tparstr], 2, PX);
                    menustr = sprintf('%s|Define %s [ %s ]', menustr, Tparstr, Tstr);                                   menuact = [ menuact 14 ];
                case 'DECTRE'
                    Lfparstr = 'Leafness parameter(s)'; [Lfstr, n_pars(end+1)] = nk_ConcatParamstr( Leafdefs );
                    PX = nk_AddParam(Leafdefs, ['ML-' Lfstr], 2, PX);
                    menustr = sprintf('%s|Define %s [ %s ]', menustr, Lfparstr, Lfstr);                                 menuact = [ menuact 16 ];
                case 'RNDFOR'
                    Dtparstr = 'Number of decision trees'; [Dtstr, n_pars(end+1)] = nk_ConcatParamstr( Treedefs );
                    PX = nk_AddParam(Leafdefs, ['ML-' Dtstr], 2, PX);
                    menustr = sprintf('%s|Define %s [ %s ]', menustr, Dtparstr, Dtstr);                                 menuact = [ menuact 17 ];
                    Dnparstr = 'Number of features'; [Dnstr, n_pars(end+1)] = nk_ConcatParamstr( NumDdefs );
                    PX = nk_AddParam(NumDdefs, ['ML-' Dnstr], 2, PX);
                    menustr = sprintf('%s|Define %s [ %s ]', menustr, Dnparstr, Dnstr);                                 menuact = [ menuact 18 ];
                case 'SEQOPT'
                    CutOffparstr = 'No. of threshold for ambiguous case propagation'; 
                    [CutOffstr, n_pars(end+1)] = nk_ConcatParamstr( SEQOPTstepsdefs );
                    PX = nk_AddParam(SEQOPTstepsdefs, ['ML-' CutOffstr], 2, PX);
                    menustr = sprintf('%s|Define %s [ %s ]', menustr, CutOffparstr, CutOffstr);                         menuact = [ menuact 22 ];
                    LimsLparstr = 'Lower population percentage(s) (- from anchor) for ambiguous case propagation'; 
                    [LimsLstr, n_pars(end+1)] = nk_ConcatParamstr( SEQOPTlimsLdefs );
                    PX = nk_AddParam(SEQOPTlimsLdefs, ['ML-' LimsLstr], 2, PX);
                    menustr = sprintf('%s|Define %s [ %s ]', menustr, LimsLparstr, LimsLstr);                           menuact = [ menuact 23 ];
                    LimsUparstr = 'Upper population percentage(s) (+ from anchor) for ambiguous case propagation'; 
                    [LimsUstr, n_pars(end+1)] = nk_ConcatParamstr( SEQOPTlimsUdefs );
                    PX = nk_AddParam(SEQOPTlimsUdefs, ['ML-' LimsUstr], 2, PX);
                    menustr = sprintf('%s|Define %s [ %s ]', menustr, LimsUparstr, LimsUstr);                           menuact = [ menuact 24 ];
                case 'WBLCOX'
                    CoxCutOffparstr = 'Percentile thresholds for class assignment';
                    [CutOffstr, n_pars(end+1)] = nk_ConcatParamstr( CoxCutoffsdefs );
                    PX = nk_AddParam(CoxCutoffsdefs, ['ML-' CutOffstr], 2, PX);
                    menustr = sprintf('%s|Define %s [ %s ]', menustr, CoxCutOffparstr, CutOffstr);                         menuact = [ menuact 25 ];  
            end
            
            GRD.GridMaxType = GridMaxType;
            if prod(n_pars) > 1 
                if OptRegul.flag, regstr = 'Yes'; else regstr = 'No'; end
                menustr = sprintf('%s|Enable regularization of model selection [ %s ]', menustr, regstr);               menuact = [ menuact 7 ];
                if OptRegul.flag
                    if isfield(TrainParam,'RFE') && ...
                            (isfield(TrainParam.RFE,'Filter') && ...
                            isfield(TrainParam.RFE.Filter,'EnsembleStrategy') && ...
                            isfield(TrainParam.RFE.Filter.EnsembleStrategy,'type') && ...
                            TrainParam.RFE.Filter.EnsembleStrategy.type ~= 9 || ...
                        isfield(TrainParam.RFE,'Wrapper') && ...
                            isfield(TrainParam.RFE.Wrapper,'EnsembleStrategy') && ...
                            isfield(TrainParam.RFE.Wrapper.EnsembleStrategy,'type') && ...
                            TrainParam.RFE.Wrapper.EnsembleStrategy.type ~= 9)
                        switch OptRegul.flag
                            case 1
                                crossstr = 'Model complexity';
                            case 2
                                crossstr = 'Ensemble diversity';
                            case 3
                                crossstr = 'Mixed Criterion (Model Complexity & Ensemble diversity)';
                        end
                        menustr = sprintf('%s|Criterion for cross-parameter model selection [ %s ]', menustr, crossstr);  menuact = [ menuact 8 ];
                    else
                        GRD.OptRegul.type = OptRegul.type;
                    end
                    menustr = sprintf('%s|Define weight (lambda) of SV ratio [ %g ]', menustr, OptRegul.lambda);           menuact = [ menuact 9 ];
                    menustr = sprintf('%s|Define non-linearity (big gamma) of SV ratio [ %g ]', menustr, OptRegul.big_gamma);  menuact = [ menuact 10 ];

                end
                switch NodeSelect.mode
                    case 1
                        modelselstr = 'Single optimum model (no ensemble)';
                    case 2
                        modelselstr = sprintf('Aggregated cross-parameter ensemble %s %g-percentile', CompStr, NodeSelect.perc);
                    case 3
                        modelselstr = sprintf('Optimized cross-parameter ensemble %s %g-percentile', CompStr, NodeSelect.perc);
                    case 4
                        modelselstr = sprintf('Variable-threshold cross-parameter ensemble');
                end
                menustr = sprintf('%s|Specify cross-parameter model selection process [ %s ]', menustr, modelselstr);     menuact = [ menuact 11 ]; 
            end
            % =============================================================
            nk_PrintLogo
            if prod(n_pars)>3, fprintf('\n');cprintf('blue','-> No. of ML parameter combinations: %g ', prod(n_pars)); end
            fprintf('\n\n'); mestr = 'Model optimization parameters'; navistr = [parentstr ' >>> ' mestr]; cprintf('*blue','You are here: %s >>> ',parentstr); 
            act = nk_input(mestr, 0,'mq', menustr, menuact );
            
            switch act
                case 1
                    Cdefs         = nk_input([Cparstr ' range'],0,'e',Cdefs);                           PX = nk_AddParam(Cdefs, ['ML-' Cparstr], 2, PX);
                case 2
                    Epsdefs       = nk_input([Rparstr ' range'],0,'e',Epsdefs);                         PX = nk_AddParam(Epsdefs, ['ML-' Rparstr], 2, PX);
                case 3
                    Nudefs        = nk_input([Rparstr ' range [0 ... 1]'],0,'e',Nudefs);                PX = nk_AddParam(Nudefs, ['ML-' Rparstr], 2, PX);
                case 4
                    Gdefs         = nk_input([Gparstr ' range'],0,'e',Gdefs);                           PX = nk_AddParam(Gdefs, ['ML-' Gparstr], 2, PX);
                case 5
                    PolyCoefdefs  = nk_input([Pcparstr ' for polynomial kernels'],0,'e',PolyCoefdefs);  PX = nk_AddParam(PolyCoefdefs, ['ML-' Pcparstr], 2, PX);
                case 6
                    PolyDegrdefs  = nk_input([Pdparstr ' for polynomial kernels'],0,'e',PolyDegrdefs);  PX = nk_AddParam(PolyDegrdefs, ['ML-' Pdparstr], 2, PX);
                case 7
                    if ~OptRegul.flag, OptRegul.flag = 2; end
                    OptRegul.flag   = nk_input('Enable model selection across CV1 parameters ?', 0, ...
                                                    'yes|no',[1,0], OptRegul.flag);
                case 8
                    OptRegul.type   = nk_input('Regularize with',0,'mq', ...
                                              ['Model complexity|' ...
                                               'Ensemble diversity|' ...
                                               'Mixed Criterion (Model Complexity & Ensemble diversity)'],1:3,OptRegul.type);
                    switch OptRegul.type 
                        case 1
                            OptRegul.RegulTypeComplexity = nk_RegulFunc_config;
                        case 2
                            OptRegul.RegulTypeDiversity = nk_RegulFunc_config;
                        case 3
                            OptRegul.RegulTypeComplexity = nk_RegulFunc_config;
                            OptRegul.RegulTypeDiversity = nk_RegulFunc_config;
                    end
                case 9
                    OptRegul.lambda = nk_input('Lambda (weight of SV ratio in score => 1 = equal weight): ', ...
                                               0, 'e', OptRegul.lambda, 1);
                case 10                
                    OptRegul.big_gamma = nk_input('Big gamma (non-linearity for SV ratio (.5 = sqrt, 1 = lin, 2 = square): ', ...
                                               0, 'e', OptRegul.big_gamma, 1);
                case 11
                    Nodemode = nk_input('Specify optimum selection mode', 0, 'mq', ...
                                               ['Select a single optimum parameter node|' ...
                                                'Generate cross-node ensemble by aggregating base learners above predefined percentile|' ... %'Generate cross-node ensemble by applying ensemble strategy on base learners above predefined threshold|' ...
                                                'Automatically determine optimal percentile for optimum cross-node ensemble performance'], [1,2,4], NodeSelect.mode);
                    if Nodemode, NodeSelect.mode = Nodemode; end
                    switch NodeSelect.mode
                        case 2
                            NodeSelect.perc = nk_input(sprintf('Models will be selected %s X%% threshold',CompStr), 0, 'e', NodeSelect.perc, 1);
                        case 3 % This option is not implemented yet in in nk_ModelNodeSelector!!!
%                             NodeSelect.perc = nk_input('Define percentage cutoff (>=%)', 0, 'e', NodeSelect.perc, 1);
%                             if ~isfield(NodeSelect,'EnsembleStrategy'), 
%                                 NodeSelect = nk_CostFun_config(NodeSelect, NM, 1);                      % Default target population
%                                 NodeSelect.SubSpaceFlag = 0;
%                                 NodeSelect = nk_SubSpaceStrategy_config(NodeSelect, NM, 1);
%                                 NodeSelect = nk_EnsembleStrategy2_config(NodeSelect, NM, 1); 
%                             end
%                             NodeSelect = nk_EnsembleStrategy2_config(NodeSelect, NM, [], navistr);
                        case 4
                            NodeSelect.percvec = nk_input('Define [lower : stepping : upper ] percentile search vector ', 0, 'e', NodeSelect.percvec);
                    end
                case 12
                    Neurondefs    = nk_input([Neuronparstr ' range'],0,'e',Neurondefs);                 PX = nk_AddParam(Neurondefs, ['ML-' Neuronparstr], 2, PX);
                case 13
                    Kdefs         = nk_input([Kparstr ' range'],0,'e',Kdefs);                           PX = nk_AddParam(Kdefs, ['ML-' Kparstr], 2, PX);
                case 14
                    Toldefs       =  nk_input([Tparstr ' range'],0,'e',Toldefs);                        PX = nk_AddParam(Toldefs, ['ML-' Tparstr], 2, PX);
                case 15
                    Weightdefs    =  nk_input([Weightparstr ' range'],0,'e',Weightdefs);                PX = nk_AddParam(Weightdefs, ['ML-' Weightparstr], 2, PX);
                case 16
                    Leafdefs    =  nk_input([Lfparstr ' range'],0,'e',Leafdefs);                        PX = nk_AddParam(Leafdefs, ['ML-' Lfparstr], 2, PX);
                case 17
                    Treedefs    =  nk_input([Dtparstr ' range'],0,'e',Treedefs);                        PX = nk_AddParam(Treedefs, ['ML-' Dtparstr], 2, PX);
                case 18
                    NumDdefs    =  nk_input([Dnparstr ' range'],0,'e',NumDdefs);                        PX = nk_AddParam(NumDdefs, ['ML-' Dnparstr], 2, PX);
                case 20
                    GRD.matLearn = nk_matLearn_config(GRD.matLearn,SVM.matLRN.learner.framework,2);
                    if isfield(GRD.matLearn,'Params')
                        for j=1:numel(GRD.matLearn.Params)
                            PX = nk_AddParam(GRD.matLearn.Params(j).range, ['ML-' GRD.matLearn.Params(j).name], 2, PX);
                        end
                    end
                case 21
                     if isfield(GRD,(SVM.prog)), PXX = GRD.(SVM.prog); else, PXX=[]; end
                    switch SVM.prog
                        case {'GLMNET','GRDBST'}
                            GRD.(SVM.prog) = nk_GLMNET_config(SVM.prog,PXX, 0);
                        case 'ROBSVM'
                            t_act = 1; while t_act > 0, [ PXX, t_act ] = nk_ROBSVM_config(SVM.prog,PXX, NM.modeflag,0); end
                            GRD.(SVM.prog) = PXX; 
                    end
                case 22
                    SEQOPTstepsdefs =  nk_input([CutOffparstr ' range'],0,'e',SEQOPTstepsdefs);         PX = nk_AddParam(SEQOPTstepsdefs, ['ML-' CutOffparstr], 2, PX);
                case 23
                    SEQOPTlimsLdefs =  nk_input([LimsLparstr ' range'],0,'e',SEQOPTlimsLdefs);          PX = nk_AddParam(SEQOPTlimsLdefs, ['ML-' LimsLparstr], 2, PX);
                case 24
                    SEQOPTlimsUdefs =  nk_input([LimsLparstr ' range'],0,'e',SEQOPTlimsUdefs);          PX = nk_AddParam(SEQOPTlimsUdefs, ['ML-' LimsUparstr], 2, PX);
                case 25
                    CoxCutoffsdefs =  nk_input([CoxCutOffparstr ' range'],0,'e',CoxCutoffsdefs);        PX = nk_AddParam(CoxCutoffsdefs, ['ML-' CoxCutOffparstr], 2, PX);
            end
            if ~isempty(PX) && ~isempty(PX.opt), n_pars = size(PX.opt,1); else n_pars = 0; end
        else
            n_pars = 0;
        end
        GRD.Cparams             = Cdefs;
        GRD.Gparams             = Gdefs;
        GRD.Epsparams           = Epsdefs;
        GRD.Nuparams            = Nudefs;
        GRD.PolyCoefparams      = PolyCoefdefs;
        GRD.PolyDegrparams      = PolyDegrdefs;
        GRD.Neuronparams        = Neurondefs;
        GRD.Kparams             = Kdefs;
        GRD.Tolparams           = Toldefs;
        GRD.Weightparams        = Weightdefs;
        GRD.Leafparams          = Leafdefs;
        GRD.Treeparams          = Treedefs;
        GRD.NumDparams          = NumDdefs;
        GRD.CutOffparams        = SEQOPTstepsdefs;
        GRD.LimsLparams         = SEQOPTlimsLdefs;
        GRD.LimsUparams         = SEQOPTlimsUdefs;
        GRD.CoxCutoffparams     = CoxCutoffsdefs;
        GRD.OptRegul            = OptRegul;
        GRD.NodeSelect          = NodeSelect;
        GRD.n_params            = n_pars;
    case 2
        
        % ***************** Setup for simulated annealing *****************
        
        % Model complexity influence:
        % ===========================
        GRD.lambda  = ...
            nk_input(['Lambda: influence of model complexity' ...
                    '(lambda=1: CV accuracy and complexity have same weigth)'],0,'e',1);
                
        GRD.big_gamma = nk_input('Penalizing gamma',0,'e',0.5);
        
        % SA params:
        % ==========
        GRD.init_T  = nk_input('Initial temperature',0,'e',100);
        GRD.end_T   = nk_input('End temperature',0,'e',1);
        GRD.dt      = nk_input('Default drop in temperature (0.1=fast; 0.01=slow)',0,'e',0.01);
        
        % Grid range params:
        % ==================        
        if isfield(GRD,'Cparams')
            defs = GRD.Cparams';
        else
            defs = [-5 2 15];
        end

        GRD.Cparams =  nk_input('C parameter exponent => St:Step:End',0,'e',defs,[3,1]);

        if isfield(GRD,'Gparams')
            defs = GRD.Gparams';
        else
            defs = [-15 2 3];
        end

        switch NM.SVM.kernel.typ
            case ' -t 1'
            case {' -t 2', ' --t 2'}
                GRD.Gparams =  nk_input('Gamma parameter exponent => St:Step:End',0,'e',defs,[3,1]);
            case ' -t 3'
        end
        
end

GRD.opt_probability = false;
GRD.PX = PX; 
TrainParam.GRD = GRD;

