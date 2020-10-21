
% ==========================================================================
function param = nk_matLearn_DefAlgoParams_config(param, options, mode, ind)

if ~exist('ind','var') || isempty(ind), ind=0; end
if exist('mode','var') && ~isempty(mode) && mode
    mode = true;
    [mn_str, mn_act] = nk_matLearn_BuildMenu_config(param, options);
    if isempty(mn_str), return; end
    nk_PrintLogo
    ind = nk_input(sprintf('Select which parameter to modify for %s',char(param.algo)),0,'mq',mn_str,mn_act);
    if ~ind, return; end
else
    mode = false;
    ind = ind + 1; 
end
[act, param] = nk_matLearn_IO_config(options, param, ind, mode);
if ~isfield(param,'Params'); param.Params=[]; end
if (ischar(act) && strcmp(act,'BACK')) || isempty(act), return; end
if strcmp(options.name{ind},'kernelOptions')
    cnt=0;
    if iscell(param.kernelOptions.name)
        for i=1:numel(param.kernelOptions.name)
            param.Params(ind+cnt).range = act.Params(i).range;
            param.Params(ind+cnt).name = sprintf('%s.%s',options.name{ind},act.Params(i).name);
            cnt=cnt+1;
        end
    else
        param.Params(ind+cnt).range = act;
        param.Params(ind+cnt).name = sprintf('%s.%s',options.name{ind},param.kernelOptions.name);
    end
    ind = ind+cnt;
elseif strcmp(options.name{ind},'subOptions')
    cnt=0;
    if iscell(param.sublearner.options.name)
        for i=1:numel(param.sublearner.options.name)
            param.Params(ind+cnt).range = act(i).range;
            param.Params(ind+cnt).name = sprintf('%s.%s',options.name{ind},act(i).name);
            cnt = cnt+1;
        end
    end
    ind = ind+cnt;
else
    param.Params(ind).range = act;
    param.Params(ind).name = options.name{ind}; 
end
%if ind <= numel(options.desc)
    param = nk_matLearn_DefAlgoParams_config(param, options, mode, ind);
%end