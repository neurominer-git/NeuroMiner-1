function [str,n_pars] = nk_ConcatParamstr(param, longflg)

str =  []; if ~exist('longflg','var') || isempty(longflg), longflg = 0; end
n_pars = numel(param);
if longflg
    for i=1:n_pars
        str = conv_param(param(i),str);
    end
else
    if n_pars > 1
        str = sprintf('%g Params: ',n_pars);
        str = sprintf('%s%s (first) -> %s (last)', str, conv_param(param(1)),conv_param(param(end)));
    elseif n_pars == 1
        str = conv_param(param);
    end
end

function str = conv_param(param, str)

flg = true;
if ~exist('str','var') || isempty(str), flg = false; end
    
if iscell(param)
    ip = param{1};
else
    ip = param;
end

if flg
    if ischar(ip)
        str = sprintf('%s|%s', str, ip);
    elseif isnumeric(ip)
        str = sprintf('%s|%g', str, ip);
    elseif isempty(ip)
        str = sprintf('%s|%s', str, 'undefined');
    else
        str = 'unknown parameter type';
    end
else
    if ischar(ip)
        str = sprintf('%s', ip);
    elseif isnumeric(ip)
        str = sprintf('%g', ip);
    elseif isempty(ip)
        str = 'undefined';
    else
        str = 'unknown parameter type';
    end
    
end