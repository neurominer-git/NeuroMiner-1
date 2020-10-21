function param = nk_Simba_config(param, setupfl)
% function param = nk_Simba_config(param, setupfl)
%
% Setup parameters for the Simba algorithm
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) N. Koutsouleris 03/2010

max_iter = 5;
block_size = 1;
utilfunc = 1;
iter_abort_crit=100;

if ~setupfl
    if isfield(param,'simba')
        max_iter = param.simba.extra_param.max_iter;
        block_size = param.simba.extra_param.block_size;
        utilfunc = param.simba.utilfunc;
        iter_abort_crit = param.simba.extra_param.iter_abort_crit;
    end
end
extra_param.max_iter = nk_input('Maximum iterations',0,'e',max_iter);
extra_param.block_size = ...
    nk_input('Block size (has to be less than subjects per CV1 training sample)',0,'e',block_size);
param.simba.utilfunc = nk_input('Utility function',0,'linear|sigmoid',[1,2],utilfunc);
if param.simba.utilfunc == 2
    betadef = nk_input('Beta of sigmoid function',0,'m','Auto-detect beta at computation|Define beta now',1:2,1);
    switch betadef
        case 1
            extra_param.beta = 'auto';
        case 2
            extra_param.beta = nk_input('Beta',0,'e',1);
    end
end
extra_param.iter_abort_crit = nk_input('Break iterations at feature stability criterion [%]', ...
    0,'e',iter_abort_crit);
param.simba.extra_param = extra_param;
param.simba.gpu = nk_input('CUDA available?',0,'yes|no',[1,0]);