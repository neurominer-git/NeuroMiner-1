function options = nk_GenMatLearnOptions(Params)

options = [];
if ~isempty(Params.desc)
    for i=1:numel(Params.val)
        if iscell(Params.val)
            options = build_option(options,Params.desc{i},Params.val{i});
        else
            options = build_option(options,Params.desc{i},Params.val(i));    
        end
    end
end

function options = build_option(options,desc,val)

if strcmp(val(1),'@'), val = str2func(val); end
str = strsplit(desc,'.');
if ~isa(val,'function_handle') 
    if ~isfinite(val), val=[]; end
end
if numel(str)>1
    options.(str{1}).(str{2}) = val; 
else
    options.(desc) = val;
end