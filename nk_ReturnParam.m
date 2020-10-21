function [Params, Desc] = nk_ReturnParam(Desc, Params_desc, opt)

Params = []; 
if iscell(Params_desc)
    for p=1:numel(Params_desc)
        switch Params_desc{p}
            case Desc
                if iscell(opt(p))
                    Params = opt{p};
                else
                    Params = opt(p);
                end
        end
    end
else
    Params = opt;
end
%if isempty(Params)
%    error(['No ' Desc ' parameters found in parameter array!!!']);
%end