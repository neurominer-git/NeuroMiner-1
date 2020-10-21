function param = nk_MKLRVM_config(param,setupfl,defaultsfl)
% =========================================================================
% FORMAT param = nk_MKLRVM_config(param, setupfl, defaultsfl)
% =========================================================================
%
% Setup parameters for RVM algorithm by Psorakis and Damoulas
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) N. Koutsouleris 03/2012

% Defaults
converg = 1;
plot_flag = 0;
standardize_flag = 1;
nmax = [];
algo = 'RVM1';
funcname_learn = 'train_mRVM1';
funcname_predict = 'predict_mRVM';
if ~exist('setupfl','var'), setupfl=0; end
if ~exist('defaultsfl','var'), defaultsfl=0; end

if ~defaultsfl
    
    if ~setupfl
        if isfield(param,'MKLRVM')
            if isfield(param.MKLRVM,'converg'), converg = param.MKLRVM.converg; end
            if isfield(param.MKLRVM,'plot_flag'), plot_flag = param.MKLRVM.plot_flag; end
            if isfield(param.MKLRVM,'standardize_flag'), standardize_flag = param.MKLRVM.standardize_flag; end
            if isfield(param.MKLRVM,'nmax'), nmax = param.MKLRVM.nmax; end
            if isfield(param.MKLRVM,'RVMalgo'), algo = param.MKLRVM.RVMalgo; end;
        end
    end
    
    if plot_flag, plot_flag_str = 'enabled'; else plot_flag_str='disabled';end
    switch converg 
        case 1
            conv_str = 'Trivial change in hyperpriors';
        case 2
            conv_str = 'Trivial change in hyperpriors & minimum number of iterations N';
        case 3
            conv_str = 'Max # of iterations';
    end
    
    if standardize_flag, standardize_flag_str = 'enabled'; else standardize_flag_str ='disabled';end
    if isempty(nmax), nmax_str = '# of training samples'; else nmax_str = num2str(nmax);end
    if strcmp(algo,'RVM1'), algonum=1;else algonum=2; end
    % ---------------------------------------------------------------------
    nk_PrintLogo
    act = nk_input('MKL-RVM algorithm: Parameter Setup',0,'mq', ...
        ['Convergence criterion [ ' conv_str ' ]|' ...
        'Plotting flag [ ' plot_flag_str ' ]|' ...
        'Standardize data [ ' standardize_flag_str ' ]|' ...
        '# of iterations [ ' nmax_str ' ]|' ...
        'RVM algorithm [ ' algo ' ]|'],1:5);
    
    switch act
        case 1
            converg = nk_input('Convergence criterion',0,'mq',...
                ['Trivial change in hyperpriors|'...
                'Trivial change in hyperpriors & minimum number of iterations N|' ...
                'Max # of iterations'],1:3,converg);
        case 2
            plot_flag = nk_input('Plot results',0,'yes|no',[1,0],plot_flag);
        case 3
            standardize_flag = nk_input('Standardize data',0,'yes|no',[1,0],standardize_flag);
        case 4
            nmaxflag = nk_input('Specify max # of iterations',0,'yes|no',[1,0]);
            if nmaxflag, 
                nmax = nk_input('Max # of iterations',0,'e',nmax);
            else
                nmax = [];
            end
        case 5
            algo = nk_input('RVM algorithm',0,'RVM1|RVM2',['RVM1';'RVM2'],algonum);
            switch param.MKLRVM.RVMalgo
                case 'RVM1'
                    funcname_learn = 'train_mRVM1';
                    funcname_predict = 'predict_mRVM';
                case 'RVM2'
                    funcname_learn = 'train_mRVM2';
                    funcname_predict = 'predict_mRVM';
            end
    end
end
param.MKLRVM.converg = converg;
param.MKLRVM.funcname_learn = funcname_learn;
param.MKLRVM.funcname_predict = funcname_predict;
param.MKLRVM.plot_flag = plot_flag;
param.MKLRVM.standardize_flag = standardize_flag;
param.MKLRVM.nmax = nmax;
param.MKLRVM.RVMalgo = algo;
if act, param = nk_MKLRVM_config(param,setupfl,defaultsfl); end

return