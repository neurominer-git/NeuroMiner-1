function ind = isnanincell(a)

ind = cellfun(@(x) any(isnan(x)),a);

