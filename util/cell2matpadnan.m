function CC = cell2matpadnan(C, pad, maxLengthCell)

if ~exist('pad','var') || isempty(pad), pad=NaN; end 
if ~exist('maxLengthCell','var') || isempty(maxLengthCell)
    maxLengthCell=max(max(cellfun('size',C,2)));  %finding the longest vector in the cell array
end
CC = cell(size(C)); n = 1; o=1;
if ndims(CC)>1, [m,n,o] = size(CC); end
for k=1:o
    for j=1:n
        for i=1:m
            CC{i,j,k} = [C{i,j,k} repmat(pad, 1, maxLengthCell - numel(C{i,j,k})) ];   %zeropad the elements in each cell array with a length shorter than the maxlength
        end
    end
end
CC=cell2mat(CC); %A is your matrix