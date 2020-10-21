function N = nk_FormatNicely(S)

[m,n] = size(S);
I = strfind(cellstr(S),'[');
ind0 = ~cellfun(@isempty,I); 

if ~isempty(I)
    try
        maxI = max(cell2mat(I(ind0))); 
    catch
        maxI = max(cellfun(@(v) v(1), I(ind0)));
    end
else
    maxI=[];
end
if isempty(maxI), N = S; return; end
N = cell(m,1);
for i=1:m
    if ~isempty(I{i})
        N{i} = sprintf(['%-' num2str( maxI ) 's%s'],S(i,1:I{i}-1),deblank(S(i,I{i}:n)));
    else
        N{i} = sprintf(['%-' num2str( maxI ) 's'],S(i,:));
    end
end
N = char(N);
