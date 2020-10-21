% ==========================================================================
function [out, param] = nk_matLearn_IO_config(opt, param, ind, mode)

out = [];
if isempty(opt) || ~isfield(opt,'desc'), return; end

if exist('ind','var') && ~isempty(ind) && iscell(opt.desc)
    if ind>numel(opt.desc), return, end
    desc = opt.desc{ind};
    def = opt.def{ind};
    sel = opt.sel{ind};
    nam = opt.name{ind};
    format = opt.format{ind};
else
    desc = opt.desc;
    def = opt.def;
    sel = opt.sel;
    nam = opt.name;
    format = opt.format;
end

if ~exist('mode','var') || isempty(mode), mode = true; end

switch format
    case 'e'
        if isfield(param,nam), def = param.(nam); end
        if mode, out = nk_input(desc,0,'e',def); else, out=def; end
    case 's'
        if isfield(param,nam), def = param.(nam); end
        if mode, out = nk_input(desc,0,'s',def); else, out=def; end
    case 'i'
        if isfield(param,nam), def = param.(nam); end
        if mode, out = nk_input(desc,0,'i',def); else, out=def; end
    case 'yes|no'
        if mode
            if isfield(param,nam),
                switch param.nam
                    case 0
                        def = 2;
                    case 1
                        def = 1;
                end
            end         
            out = nk_input(desc,0,format,sel,def);
        else
            out = def;
        end
    case 'model_selector'
        if ~isfield(param,'sublearner'), param.sublearner = []; end
        param.sublearner.framework = param.learner.framework;
        param.sublearner = nk_matLearn_getopts_config(param.sublearner,'get_sublearners',[],[],param.sublearner.framework);
        param.sublearner.algo = nk_matLearn_IO_config(param.sublearner, param, ind, mode);
        out = sprintf('@ml_%s_%s', param.sublearner.framework, char(param.sublearner.algo));
        if isfield(param,'Params')
            ind = false(1,numel(param.Params));
            if ~isempty(param.Params)
                for j=1:numel(param.Params)
                   if strfind(param.Params(j).name,'subOptions')
                       ind(j)=true;
                   end 
                end
                param.Params(ind)=[];
            end
        end
        
    case 'options_selector'
        if isfield(param,'sublearner')
            if ~isfield(param.sublearner,'options'), param.sublearner.options = []; end
            param.sublearner.options = nk_matLearn_getopts_config(param.sublearner.options,'get_learner_params', char(param.sublearner.algo),param.sublearner.framework);
            param.sublearner = nk_matLearn_DefAlgoParams_config(param.sublearner, param.sublearner.options, mode, ind);
        end
        out = param.sublearner.Params;
        
    case 'kernel_func_selector'
        if ~isfield(param,'kernelFunc'), param.kernelFunc = []; end
        param.kernelFunc = nk_matLearn_getopts_config(param.kernelFunc,'get_kernel_func',char(param.algo),param.learner.framework);
        out = nk_matLearn_IO_config(param.kernelFunc, param, ind, mode);
        ind = false(1,numel(param.Params));
        if ~isempty(param.Params)
            for j=1:numel(param.Params)
               if strfind(param.Params(j).name,'kernelOptions')
                   ind(j)=true;
               end 
            end
            param.Params(ind)=[];
        end
        
    case 'kernel_params_selector'
        if ~isfield(param,'kernelOptions'), param.kernelOptions = []; end
        indx = 0;
        for i=1:numel(param.Params)
            if strcmp(param.Params(i).name,'kernelFunc'); indx = i; break; end
        end
        if indx, 
            param.kernelOptions = nk_matLearn_getopts_config(param.kernelOptions,'get_kernel_params',char(param.Params(indx).range),param.learner.framework);
            if isfield(param.kernelOptions,'name') && ~ischar(param.kernelOptions.name)
                out = nk_matLearn_DefAlgoParams_config(param.kernelOptions, param.kernelOptions, mode);
            else
                out = nk_matLearn_IO_config(param.kernelOptions, param, indx, mode);
            end
        else
            out=[];
        end
        
    otherwise
        if isfield(opt,nam), def = find(strcmp(sel,opt.(nam))); else def=1;end
        if mode, out = nk_input(desc,0,'mq',format,sel,def); else, out = sel{def}; end
end