function [Dummy, Fu] = nk_MakeDummyVariables(V,mode)

if ~exist('mode','var'), mode = 'stable'; end
[m, n] = size(V);
if m > 1 && n > 1, error('Only vector operations are supported!'); end
[Fu, IA, CA] = unique(V,mode); nFu = numel(Fu);
if m > 1
    Dummy = false(m, nFu);
else
    Dummy = false(nFu, n);
end

for i = 1:nFu
    try
        ind = V == Fu(i);
    catch
        ind = strcmp(V,Fu{i});
    end
    if m > 1 
        Dummy(ind,i) = true;
    else
        Dummy(i, ind) = true;
    end
end
