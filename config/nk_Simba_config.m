function param = nk_Simba_config(param, setupfl, defaultsfl, ngroups)
% function param = nk_Simba_config(param, setupfl)
%
% Setup parameters for Simba algorithm
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) N. Koutsouleris 04/2015

% Defaults
max_iter = 5;
block_size = 1;
st_points = 1;
utilfunc = 1;
utilfuncstr = 'linear'; 
Beta = 'auto'; betadef = 1; betastr='';
cuda = 0;
iter_abort_crit=100;

if ngroups > 2, binmode = 0; else binmode = 1; end

if ~exist('setupfl','var') || isempty(setupfl), setupfl=0; end
if ~exist('defaultsfl','var') || isempty(defaultsfl), defaultsfl=0; end

if ~defaultsfl
    if ~setupfl
        if isfield(param,'simba')
            if isfield(param.simba,'extra_param')
                if isfield(param.simba.extra_param,'max_iter'), ...
                        max_iter = param.simba.extra_param.max_iter; end
                if isfield(param.simba.extra_param,'block_size'), ...
                        block_size = param.simba.extra_param.block_size; end
                if isfield(param.simba.extra_param,'start_points'), ...
                        st_points = param.simba.extra_param.start_points; end
                if isfield(param.simba.extra_param,'iter_abort_crit'), ...
                        iter_abort_crit = param.simba.extra_param.iter_abort_crit; end
                if isfield(param.simba.extra_param,'beta')...
                        Beta = param.simba.extra_param.beta; end
            end
            if isfield(param.simba,'utilfunc')
                utilfunc = param.simba.utilfunc;
                if utilfunc == 2
                    if isfield(param.simba.extra_param,'beta'), 
                        betastr = param.simba.extra_param.beta;                    
                    else
                        betastr = Beta;
                    end
                    if isnumeric(betastr), 
                        Beta = betastr;
                        betastr = nk_ConcatParamstr(betastr);
                        betadef = 2;
                    else
                        betadef = 1;
                    end
                    betastr = [' (beta: ' betastr ')'];
                else
                    betastr = '';
                end
                utilfuncstr = param.simba.extra_param.utility;
            end
            if isfield(param.simba,'cuda'), cuda = param.simba.gpu; end
        end
    end

    %if ~cuda, cudastr = 'no'; else cudastr = 'yes'; end
    if isfield(param,'binmode'), binmode = param.binmode; end
    if ~binmode, binmodestr='multi-group'; else binmodestr='pairwise'; end

    % -------------------------------------------------------------------------
    nk_PrintLogo
    if ngroups > 2,

        act = nk_input('Simba: Parameter Setup',0, 'mq', ...
            ['# Starting points to avoid local maxima [' num2str(st_points) ']|' ...
            'Maximum # of iterations [' num2str(max_iter) ']|' ...
            'Iteration stopping criterion [' num2str(iter_abort_crit) '% feature rank stability]|' ...
            'Block size [' num2str(block_size) ']|' ...
            'Utility function [' utilfuncstr betastr ']|' ...
            'Multi-group / Pairwise group processing [' binmodestr ']'],1:6);
    else
        act = nk_input('Simba: Parameter Setup',0, 'mq', ...
            ['# Starting points to avoid local maxima [' num2str(st_points) ']|' ...
            'Maximum # of iterations [' num2str(max_iter) ']|' ...
            'Iteration stopping criterion [' num2str(iter_abort_crit) ']|' ...
            'Block size [' num2str(block_size) ']|' ...
            'Utility function [' utilfuncstr betastr ']'],1:5);
    end

    switch act
        case 1
            st_points = nk_input('Starting points',0,'i',st_points,1);
        case 2
            max_iter = nk_input('Maximum iterations',0,'i',max_iter,1);
        case 3
            iter_abort_crit = nk_input('Break iterations at feature stability criterion [%]', ...
                0,'e',iter_abort_crit,1);
        case 4
            block_size = ...
                nk_input('Block size (has to be less than subjects per CV1 training sample)',0,'i', ...
                block_size,1);
        case 5
            utilfunc = nk_input('Utility function',0,'linear|sigmoid', ...
                [1,2],utilfunc);

            switch utilfunc
                case 1
                    utilfuncstr = 'linear';
                    Beta = 'auto';
                case 2
                    utilfuncstr = 'sigmoid';
                    betadef = nk_input('Beta of sigmoid function',0,'m', ...
                        'Auto-detect beta at computation|Define beta now',1:2, betadef);
                    switch betadef
                        case 1
                            Beta = 'auto';
                        case 2
                            Beta = nk_input('Beta',0,'e', Beta);
                    end
            end
        case 6
            binmode = uint8(nk_input('Multi-group or pairwise processing?', 0, 'm', ...
                                'Multi-group|Binary',[0,1],binmode));
    end
else
    act = 0;
end
param.simba.extra_param.start_points = st_points;
param.simba.extra_param.max_iter = max_iter;
param.simba.extra_param.iter_abort_crit = iter_abort_crit;
param.simba.extra_param.block_size = block_size;
param.binmode = binmode;
param.simba.gpu = cuda;
param.simba.utilfunc = utilfunc;
param.simba.extra_param.beta = Beta;
param.simba.extra_param.utility = utilfuncstr;

if act, param = nk_Simba_config(param, setupfl, defaultsfl, ngroups); end
