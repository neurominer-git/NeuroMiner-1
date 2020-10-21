function indx=findrow(A,a)
% function indx=findrow(A,a)
% Finds the rows in A that are a 
% returns row-indices 
indx=[];
[numrow,numcol]=size(A);
[vec,veccol]=size(a);
if (isempty(A))
    return;
end;
if (vec~=1 | veccol~=numcol)
    error('a has to be a row of A');
end;
indx=find(all(A==repmat(a,numrow,1),2));