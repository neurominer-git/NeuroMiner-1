function [I,M]=split_data(X,num_quant,varargin)
% functions [I,M]=split_data(X,numquant):
% splits the data into x quantiles and
% OUTPUT:
%   I: number of quantiles for each data point
%   M: center of each bin for each observation 
% VARARGIN:
% 'subset': only use a certain subset
% 'split': calculate the quantiles after splitting the data into different
%           categories
% v1.0 01/06/06: subset option introduced
%  1.1 18/3/08: split version introduced
subset=[];
split=[];
vararginoptions(varargin,{'subset','split'});
if (isempty(subset));
    subset=ones(size(X,1),1); 
end;
if (isempty(split))
    split=ones(size(X,1),1);
end;

cat=unique(split(find(subset),:),'rows');

% Define the quantiles 
quant=[0:num_quant-1]./num_quant;
I=X*NaN; 
M=X*NaN; 
for i=1:size(cat,1);
    indx=findrow([split subset],[cat(i,:) 1]);
    borders=prctile(X(indx,:),quant*100);
    for j=1:num_quant
        k=find(X(indx)>=borders(j));
        I(indx(k),1)=j;
    end;
    for j=1:num_quant
        k=find(I(indx)==j);
        M(indx(k),1)=mean(X(indx(k)));
    end;
end;
