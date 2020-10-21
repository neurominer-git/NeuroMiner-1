function T=field_matrix2vector(T,fieldname);
% Makes a matrix in a structure to n vectors 
% with appended numbers 
%
% function T=field_matrix2vector(T,fieldname);
% 
NAMES=fieldnames(T);
thisfield=strmatch(fieldname,NAMES,'exact');
if (isempty(thisfield))
    error('field does not exist');
end;

M=getfield(T,fieldname);
T=rmfield(T,fieldname);
% a=T.(fieldname); 
numcols=size(M,2);
if (numcols<10) 
    format='%d';
elseif (numcols<100)
    format='%2.2d';
else 
    error ('Matrices are not allowed to have more than 99 columns');
end;
for c=1:numcols 
    f=sprintf([fieldname format],c);
    if (isfield(T,f))
        error('Can not overwrite fields');
    end;
    T=setfield(T,f,M(:,c));
end;
try
    T=orderfields(T,[1:thisfield-1 length(NAMES):length(NAMES)+numcols-1 thisfield:length(NAMES)-1]);
end;