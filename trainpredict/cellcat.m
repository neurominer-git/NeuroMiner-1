function catarray = cellcat(arr, dim)

[ix,jx] = size(arr);
if nargin < 2, dim = 1; end
catarray = [];
for i=1:ix
    for j=1:jx
        switch dim
            case 1
                catarray =[catarray; arr{i,j}];
            case 2
                catarray =[catarray arr{i,j}];
        end
    end
end
