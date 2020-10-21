function param = nk_Gflip_config(param, setupfl, defaultsfl, ngroups)
% function param = nk_Gflip_config(param, setupfl)
%
% Setup parameters for G-flip algorithm
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) N. Koutsouleris 4/2015

% Defaults
max_iter = 5;
block_size = 1;
st_points = 1;
utilfunc = 2;
utilfuncstr = 'linear'; betastr='';
cuda = 0;
Beta = 'auto'; betadef = 1;
if ngroups > 2, binmode = 0; else binmode = 1; end

if ~exist('setupfl','var') || isempty(setupfl), setupfl=0; end
if ~exist('defaultsfl','var') || isempty(defaultsfl), defaultsfl=0; end

if ~defaultsfl
    
    if isfield(param,'gflip')
        if isfield(param.gflip,'extra_param')
            if isfield(param.gflip.extra_param,'max_iter'), ...
                    max_iter = param.gflip.extra_param.max_iter; end
            if isfield(param.gflip.extra_param,'block_size'), ...
                    block_size = param.gflip.extra_param.block_size; end
            if isfield(param.gflip.extra_param,'start_points'), ...
                    st_points = param.gflip.extra_param.start_points; end
        end

        if isfield(param.gflip,'utilfunc')
            utilfunc = param.gflip.utilfunc;
            if utilfunc == 3
                if isfield(param.gflip.extra_param,'beta'), 
                    betastr = param.gflip.extra_param.beta;                    
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
            utilfuncstr = param.gflip.extra_param.utility;
        end
        if isfield(param.gflip,'cuda'), cuda = param.gflip.gpu; end
    end

    %if ~cuda, cudastr = 'no'; else cudastr = 'yes'; end
    if isfield(param,'binmode'), binmode = param.binmode; end
    if ~binmode, binmodestr='multi-group'; else binmodestr='pairwise'; end

    % -------------------------------------------------------------------------
    nk_PrintLogo
    if ngroups > 2,

        act = nk_input('Gflip: Parameter Setup',0, 'mq', ...
            ['# Starting points to avoid local maxima [' num2str(st_points) ']|' ...
            'Maximum # of iterations [' num2str(max_iter) ']|' ...
            'Block size [' num2str(block_size) ']|' ...
            'Utility function [' utilfuncstr betastr ']|' ...
            'Multi-group / Pairwise group processing [' binmodestr ']'],1:5);
    else
        act = nk_input('Gflip: Parameter Setup',0, 'mq', ...
            ['# Starting points to avoid local maxima [' num2str(st_points) ']|' ...
            'Maximum # of iterations [' num2str(max_iter) ']|' ...
            'Block size [' num2str(block_size) ']|' ...
            'Utility function [' utilfuncstr betastr ']'],1:4);
    end

    switch act
        case 1
            st_points = nk_input('Starting points',0,'i',st_points,1);
        case 2
            max_iter = nk_input('Maximum iterations',0,'i',max_iter,1);
        case 3
            block_size = ...
                nk_input('Block size (has to be less than subjects per CV1 training sample)',0,'i', ...
                block_size,1);
        case 4
            utilfunc = nk_input('Utility function',0,'m','zero-one|linear|sigmoid', ...
                1:3,utilfunc);

            switch utilfunc
                case 1
                    utilfuncstr = 'zero-one';
                    Beta = 'auto';
                case 2
                    utilfuncstr = 'linear';
                    Beta = 'auto';
                case 3
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
        case 5
            binmode = uint8(nk_input('Multi-group or pairwise processing?', 0, 'm', ...
                                'Multi-group|Binary',[0,1],binmode));
    end
else
    act = 0;
end
param.gflip.extra_param.start_points = st_points;
param.gflip.extra_param.max_iter = max_iter;
param.gflip.extra_param.block_size = block_size;
param.binmode = binmode;
param.gflip.gpu = cuda;
param.gflip.utilfunc = utilfunc;
param.gflip.extra_param.beta = Beta;
param.gflip.extra_param.utility = utilfuncstr;

if act, param = nk_Gflip_config(param, setupfl, defaultsfl, ngroups);end