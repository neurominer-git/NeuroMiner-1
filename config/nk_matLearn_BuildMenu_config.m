% ==========================================================================
function [mn_str, mn_act] = nk_matLearn_BuildMenu_config(param, opt)
if isempty(opt), mn_str = []; mn_act = []; return; end
nO = numel(opt.desc);
desc = []; icnt=1;
for i=1:nO
    opt_str = '?'; no_mnu = false;
    if isfield(param,'Params') && i<=numel(param.Params) && ~isempty(param.Params(i).name)
        Pname = strsplit(param.Params(i).name,'.');
        switch Pname{1}
            case 'subModel'
                if iscell(param.sublearner.algo), algo = char(param.sublearner.algo); else algo = param.sublearner.algo; end
                opt_str = sprintf('%s', algo);
            case {'subOptions','kernelOptions'}
                jopt_str=[];
                for j=i:numel(param.Params)
                    if strfind(param.Params(j).name,Pname{1})
                        px = strsplit(param.Params(j).name,'.');
                        jopt_str = sprintf('%s%s: %s ; ', jopt_str, px{2}, nk_ConcatParamstr(param.Params(j).range));
                    end
                end
                if ~isempty(jopt_str), opt_str = jopt_str(1:end-3); end
            otherwise
                if ~isempty(opt.format{i}) && ~no_mnu
                    switch opt.format{i}
                        case 'yes|no'
                            if ~param.Params(i).range, opt_str = 'no';  else, opt_str = 'yes'; end
                        case 'e'
                            opt_str=nk_ConcatParamstr(param.Params(i).range);
                        otherwise
                            if iscell(param.Params(i).range)
                                opt_str=char(param.Params(i).range);
                            else
                                opt_str=param.Params(i).range;
                            end
                    end
                    if ~ischar(opt_str), opt_str = num2str(opt_str); end
                end
        end
        
    else
        switch opt.name{i}
            case 'subOptions'
                if isfield(param,'sublearner')
                    if iscell(param.sublearner.algo), algo = char(param.sublearner.algo); else algo = param.sublearner.algo; end
                    subopts = nk_matLearn_getopts_config([],'get_learner_params',algo);
                    if isempty(subopts), no_mnu=true; end
                end
            case 'kernelOptions'
                if isfield(param,'kernelFunc') && isfield(param,'Params')
                    for j=1:numel(param.Params)
                        if strcmp(param.Params(j).name,'kernelFunc'),
                            kernelopts = nk_matLearn_getopts_config(opt, 'get_kernel_params', char(param.Params(j).range), param);
                            if ~isfield(kernelopts,'desc'), no_mnu=true; end
                            break
                        end
                    end
                end
        end
    end
    if ~no_mnu, 
        desc{icnt} = sprintf('%s [ %s ]', opt.desc{i}, opt_str); 
        icnt=icnt+1;
    end
end
        
mn_str = strjoin(desc,'|');
mn_act = 1:numel(desc);