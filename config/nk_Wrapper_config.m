function [act, param] = nk_Wrapper_config(act, param, SVM, MODEFL, GRD, MULTI, defaultsfl, parentstr)
%
% Defaults
% --------
% Activate feature filtering
Wrapper.flag                = 0;    % Perform wrapper-based feature selection
Wrapper.optflag             = 2;    % Perform wrapper-based feature selection at parameter optimum
Wrapper.type                = 1;
Wrapper.datamode            = 2;
Wrapper.CostFun             = 2;   % currently static
[ ~, Wrapper ]              = nk_GreedySearch_config(Wrapper,SVM, MULTI, 1);
Wrapper                     = nk_SA_config(Wrapper,1); % Define SA defaults
% Subspace strategy
Wrapper.SubSpaceFlag        = 0;
Wrapper.SubSpaceStepping    = 0;
Wrapper                     = nk_EnsembleStrategy2_config(Wrapper, SVM, MODEFL, 1);
% Probabilistic Feature Extraction (PFE) across CV1 partitions 
Wrapper.PFE                 = nk_ProbalisticFea_config([],1);
if ~exist('defaultsfl','var') || isempty(defaultsfl),  defaultsfl = 0; end;

if ~isempty(act) || ~defaultsfl
    
    if isfield(param,'Wrapper'), 
        if ~isfield(param.Wrapper,'EnsembleStrategy'), 
            param.Wrapper = nk_EnsembleStrategy2_config(Wrapper, SVM, MODEFL, 1);
        end
        Wrapper = param.Wrapper; 
    end
    
    nk_PrintLogo
    mestr = 'Wrapper-based model selection setup'; navistr = [parentstr ' >>> ' mestr]; cprintf('*blue','\nYou are here: %s >>>',parentstr);
    
    if ~Wrapper.flag
        
        act = nk_input(mestr,0,'mq','Activate wrapper methods', 1, 1);
        
    else
       
        d = nk_GetParamDescription2([],param,'FeatWrap');
        
        if isfield(param,'Filter') && param.Filter.flag && isfield(param.Filter,'SubSpaceFlag') && param.Filter.SubSpaceFlag
            e = nk_GetParamDescription2([],param.Wrapper,'EnsType');
            cprintf('red','\n================================================================')
            cprintf('red','\n* Pre-Wrapper filtering/subspace optimization method detected. *')
            cprintf('red','\n* Wrapper will be applied to pre-defined feature subspaces.    *')
            cprintf('red','\n* Further subspace optimization methods enabled.               *')
            cprintf('red','\n================================================================')

            SubStr = ['Evaluate features within subspaces [ ' e.EnsStrat ' ]|'];
            actind = [1 2];    
        else
            SubStr ='';
            actind = 1;
        end
        
        if prod(GRD.n_params)>1
            GrdOptStr = ['Learn wrapper model only at parameter optimum [ ' d.WrapperOptFlag ' ]|']; actind = [actind 3 ];
        else
            GrdOptStr = '';
        end
        
        
        act = nk_input(mestr,0,'mq',['Deactivate wrapper methods|'...
                                   SubStr ...
                                   GrdOptStr ...
                                   'Wrapper type [ ' d.WrapperStr ' ]|'...
                                   'CV1 data partitions for optimization [ ' d.WrapperDataMode ' ]|'...
                                   'Cross-CV1 Feature selection [ ' d.WrapperPFE ' ]'], [actind 4:6], 1);
            
    end
    
    switch act
        
        case 1
            Wrapper.flag = ~Wrapper.flag;    
            if ~Wrapper.flag && Wrapper.optflag == 1, Wrapper.optflag = 2; end
 
        case 2
            if isfield(param,'Filter') && isfield(param.Filter,'SubSpaceFlag')
                Wrapper.SubSpaceFlag = param.Filter.SubSpaceFlag;
            else
                Wrapper.SubSpaceFlag = 0;
            end
            Wrapper = nk_EnsembleStrategy2_config(Wrapper, SVM, MODEFL, [], navistr);
            if Wrapper.SubSpaceFlag
                SubSpaceSteppingFlag = nk_input('Define feature subspace stepping',0,'yes|no',[1,0],0);
                if SubSpaceSteppingFlag
                    Wrapper.SubSpaceStepping = nk_input('Divide maximum number of feature subspaces by',0, ...
                        'w1',Wrapper.SubSpaceStepping);
                end
            end
        case 3
            if isfield(Wrapper,'optflag'), optflag = Wrapper.optflag; else, optflag = 1; end 
            Wrapper.optflag = nk_input('Learn Wrapper only at optimum',0,'yes|no',[1,2], optflag);
            
        case 4
            Wrapper.type = nk_input('Feature search mode',0,'m','Greedy feature selection|Simulated Annealing',1:2, Wrapper.type);
            switch Wrapper.type
                case 1
                    t_act = 1; while t_act>0, [ t_act , Wrapper ] = nk_GreedySearch_config(Wrapper, SVM, MULTI, 0, navistr); end
                case 2
                    Wrapper = nk_SA_config(Wrapper); 
            end
        case 5
            Wrapper.datamode = nk_input('Samples for optimization',0,'m','CV1 training data|CV1 test data|CV1 training & test data',[1,2,3], Wrapper.datamode);
        case 6
            sact = 1; while sact>0, [Wrapper.PFE, sact] = nk_ProbalisticFea_config(Wrapper.PFE, navistr); end
    end
    
end
param.Wrapper = Wrapper;
