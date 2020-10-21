function str = nk_RemCommaFromStr(str)

try
    if iscell(str)
        for n=1:numel(str)
            str{n} = regexprep(str{n},',1','');
        end
    else
        if size(str,1) > 1
            for n=1:size(str,1)
                str(n,:) = regexprep(str(n,:),',1',' ');
            end
        else
            str = regexprep(str,',1','');
        end
    end
catch
     if iscell(str)
        for n=1:numel(str)
            i = strfind(str{n},',1'); str{n}(i:end)=[];
        end
    else
        if size(str,1) > 1
            for n=1:size(str,1)
                i = strfind(str(n,:),',1');
                str(n,i:end)=[];
            end
        else
            i = strfind(str,',1');
            str(i:end)=[];
        end
    end
   
end