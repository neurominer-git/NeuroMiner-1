function param = nk_RGS_config(param, defaultsfl)
% function param = nk_RGS_config(param, setupfl)
%
% Setup parameters for the RGS algorithm (kNN-regression)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) N. Koutsouleris 04/2015

epochs = 1;
num_starts = 2;
verbose = 1;

if ~exist('defaultsfl','var'),defaultsfl=0; end

if ~defaultsfl % use interactive setup
   
    if isfield(param,'RGS')
        if isfield(param.RGS.extra_param,'epochs'),     epochs = param.RGS.extra_param.epochs; end
        if isfield(param.RGS.extra_param,'num_starts'), num_starts = param.RGS.extra_param.num_starts; end
        if isfield(param.RGS.extra_param,'verbos'),     verbose = param.RGS.extra_param.verbose; end
        if isfield(param.RGS.extra_param,'beta'),       beta = param.RGS.extra_param.beta; end
        if isfield(param.RGS.extra_param,'k'),          k = param.RGS.extra_param.k; end
    end
    
    if verbose , verbosestr = 'yes'; else verbosestr = 'no'; end
    if exist('beta','var') , 
        betastr = nk_ConcatParamstr(beta); 
        betadef = 2; 
    else
        betastr = 'automatic'; beta=1; 
        betadef = 1;
    end
    if exist('k','var') , 
        knnstr = nk_ConcatParamstr(k); 
        kdef = 2;
    else
        knnstr = 'automatic'; 
        k = 8; 
        kdef = 1;
    end

    % -------------------------------------------------------------------------
    nk_PrintLogo

    act = nk_input('RGS: Parameter Setup',0, 'mq', ...
            ['# Starting points [' num2str(num_starts) ']|' ...
            '# Epochs [' num2str(epochs) ']|' ...
            'Beta of gaussian function [' betastr ']|' ...
            '# Nearest neighbours [' knnstr ']|' ...
            'Verbose [' verbosestr ']'],1:5);

    switch act
        case 1
            num_starts = nk_input('Starting points (to avoid local maxima)',0,'i', num_starts,1);
        case 2
            epochs = nk_input('Number of epochs',0,'i',epochs,1);
        case 3
            betadef = nk_input('Beta of Gaussian function',0,'m','Auto-detect beta at computation|Define beta now',1:2,1);
            if betadef == 2
                    beta = nk_input('Beta',0,'e',beta);
            end
        case 4
            kdef = nk_input('K nearest neighbours',0,'m','Auto-detect K at computation|Define K now',1:2,1);
            if kdef == 2
                k = nk_input('# of nearest neighbours',0,'i',k);
            end
        case 5
            verbose = nk_input('Verbose',0,'yes|no',[1,0],verbose);
    end
else
    act = 0;
end
param.RGS.extra_param.num_starts = num_starts;
param.RGS.extra_param.epochs = epochs;
if exist('betadef','var') && betadef == 2
    param.RGS.extra_param.beta = beta;
else
    if isfield(param.RGS.extra_param,'beta'), ...
            param.RGS.extra_param = rmfield(param.RGS.extra_param,'beta'); end
end
if exist('kdef','var') && kdef == 2
    param.RGS.extra_param.k = k;
else
    if isfield(param.RGS.extra_param,'k'), ...
            param.RGS.extra_param = rmfield(param.RGS.extra_param,'k'); end
end
param.RGS.extra_param.verbose = verbose;

if act, param = nk_RGS_config(param, defaultsfl); end