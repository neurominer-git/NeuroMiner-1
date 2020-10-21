function D=insertrow(D,numrow,ROW)
% function D=insertrow(D,numrow,ROW)
% inserts the structure ROW at the position numrow into 
% the data-structure D
% if a field of D is not given in ROW, is is filled with NAN's 
% Joern Diedrichsen 
% Version 1.0 9/18/03
% see also setrow, getrow
field=fieldnames(D);
for f=1:length(field)
    if (isfield(ROW,field{f}))
        F=getfield(D,field{f});
        I=getfield(ROW,field{f});    
    else 
        F=getfield(D,field{f});
        I=NaN;    
    end;   
    D=setfield(D,field{f},[F(1:numrow-1);I;F(numrow:end)]);
end;
