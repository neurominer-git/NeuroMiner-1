function D=dsort(D,sort); 
% function D=dsort(D,sort); 
% Sorts a data structure in ascending order 
% of the data fields given in sort 
[Y,I]=sortrows(sort); 
D=getrow(D,I); 
