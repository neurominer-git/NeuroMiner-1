function opt_def = nk_matLearn_FindDefInParam_config(opt_format, opt_name, opt_def, opt_sel, param)
    
    switch opt_format
        case 'yes|no'
            if strcmp(param.name, opt_name) && ~param.range, opt_def = 2; else, opt_def = 1; end 
        case {'e','s','i'}
            if strcmp(param.name, opt_name), opt_def = param.range; end
        case 'model_selector'
        case 'options_selector'
        case 'kernel_func_selector'
        case 'kernel_params_selector'
        otherwise
            opt_def = find(strcmp(opt_sel,char(param.range)));
    end
    
end