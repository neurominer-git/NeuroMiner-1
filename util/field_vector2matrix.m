function T=field_vector2matrix(T,fieldname);
% function T=field_vector2matrix(T,fieldname);
% Makes a matrix in a structure to n vectors 
% with appended numbers 
NAMES=fieldnames(T);

W=[];
i=1;
while i<100;
    f=strmatch(sprintf([fieldname '%d'],i),NAMES);
    if isempty(f)
        f=strmatch(sprintf([fieldname '%2.2d'],i),NAMES);
        if isempty(f)
            break;
        end;
    end;
    fields(i)=f;
    W=[W getfield(T,NAMES{fields(i)})];
    i=i+1;
end;
if (isempty(W))
    error('field does not exist');
end;

for j=1:length(fields)
    T=rmfield(T,NAMES{fields(j)});
end;
T=setfield(T,fieldname,W);
