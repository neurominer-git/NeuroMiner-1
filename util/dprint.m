function dprint(D)
%  function dprint(filename,CA)
%       writes out a struct as a tab-delimted sheet
%       filename: if empty, function will promt user for filename
%       D: data-struct
%           each field is a column
%           function delas with empty vecotrs and ragged arrays (vectors of different length), too
%  Joern Diedrichsen
%  Department for Psychology
%  April 2002
names=fieldnames(D);
linefeed=sprintf('\n');
% put in variable names
numvar=length(names);
for v=1:numvar
    if (~strcmp(names{v},'header'))
        fprintf('%5s',names{v});
        if(v<numvar)
            fprintf('\t',linefeed);
        end;
        var_length(v)=length(getfield(D,names{v}));    
        
    end;
end;
fprintf('%s',linefeed);
% now get out all the variables 
for i=1:max(var_length)
    for v=1:numvar
        if (~strcmp(names{v},'header'))
            var=getfield(D,names{v});
            if var_length(v)>=i
                if (iscell(var))
                    fprintf('%s',var{i,1});
                elseif (ischar(var(i)))
                    fprintf('%s',deblank(var(i)));
                else 
                    fprintf('%2.2f',var(i));
                end;
            else
                %fprintf(fid,' ');          
            end;
        end;
        if v==numvar
            fprintf('%s',linefeed);
        elseif (~strcmp(names{v},'header'))
            fprintf('\t');
        end;
    end;
end;
