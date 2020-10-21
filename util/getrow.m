function ROW=getrow(D,numrow)
% function ROW=get(D,numrow)
% extracts the structure ROW at the position numrow from 
% the data-structure D
% Joern Diedrichsen 
% Version 1.0 9/18/03
% SEE also: insertrow,setrow
field=fieldnames(D);
ROW=[];
for f=1:length(field)
    F=getfield(D,field{f});
    if iscell(F)
        ROW=setfield(ROW,field{f},F(numrow,:));
    else

        ROW=setfield(ROW,field{f},F(numrow,:));
    end;
end;
