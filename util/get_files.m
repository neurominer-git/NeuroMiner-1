function P=get_files(filename)
% Resolves the * file naming convention
% returns a cell array of file names, that match the filename
% needs to be in the working directory 
indx=find(filename=='*');P={};
if (isempty(indx))
    P{1}=filename;
    return;
elseif (length(indx)>1)
    error('get_files: Only one placeholder allowed at a time');
else 
    before=indx-1;
    after=length(filename)-indx;
    LIST=dir;
    for l=1:length(LIST);
        name=LIST(l).name;
        if (length(name)>=length(filename))
            if (before==0 | strcmp(name(1:before),filename(1:before)))
                if (after==0 | strcmp(name(end-after+1:end),filename(indx+1:end)))
                    P{end+1}=[pwd '/' name]; 
                end;
            end;
        end;
    end;
end;
