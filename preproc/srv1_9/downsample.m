function D=downsample(tensor,p)
% downsample images
% tensor: 3-way tensor, tensor(:,:,i) is the ith image
% p, integer, 1/p data will be returned
% D: 3-way tensor, the downsampled data
% can be extended to n-way
% Yifeng Li
% September 07, 2011

st=size(tensor);
rows=1:p:st(1);
cols=1:p:st(2);
D=tensor(rows,cols,:);
end