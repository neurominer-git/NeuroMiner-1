function param = nk_GLMNET_config(prog, param, defaultsfl)

if ~exist('defaultsfl','var') || isempty(defaultsfl), defaultsfl = true; end
param.options = nk_matLearn_getopts_config([], 'get_learner_params', prog, param);
ind = 1;
if ~defaultsfl
    while ~isempty(ind)
        [mn_str, mn_act] = nk_matLearn_BuildMenu_config(param, param.options);
        nk_PrintLogo; 
        ind = nk_input('Select parameter to modify',0,'mq',mn_str,mn_act);
        if ~ind, break; end
        [out, param] = nk_matLearn_IO_config(param.options, param, ind);
        param.Params(ind).range = out;
        param.Params(ind).name = param.options.name{ind}; 
    end
else
    for ind = 1:numel(param.options.name)
        param.Params(ind).range = param.options.def{ind};
        param.Params(ind).name = param.options.name{ind}; 
    end
end

