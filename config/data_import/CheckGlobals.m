function [sz_str, mess, IO] = CheckGlobals(IO, act, na_str, mess)

if ~exist('mess','var') , mess=[]; end
sz_str = na_str;
if strcmp(IO.globvar_edit,na_str), return; end
nP = size(IO.PP,1);
flg=false;
try
    switch act
        case 'check_globvar'
            
            globvar_str = sprintf('Var: ''%s''',IO.globvar_edit);
            
            g = evalin('base', IO.globvar_edit);
            
            if ~isnumeric(g)
                mess = GenerateMessageEntry(mess,  sprintf('ERROR: The globals variable ''%s'' does not contain numeric data.',IO.globvar_edit));
                flg = true;
            end
            if size(g,2) > 1
                mess = GenerateMessageEntry(mess, sprintf('ERROR: The globals variable ''%s'' should be a vector, not a matrix.',IO.globvar_edit, nP));
                flg = true;
            end
            if size(g,1) ~= nP 
                mess = GenerateMessageEntry(mess, sprintf('ERROR: The globals variable ''%s'' should be a vector with exactly %g entries in the order of the images.',IO.globvar_edit, nP));
                flg = true;
            end
            ind = 1:nP;
            
        case 'check_globfile'
           
            [~,n,e] = fileparts(IO.globvar_edit);
            
            globvar_str = sprintf('File: ''%s%s''',n,e);
            GG = readtable(IO.globvar_edit,'ReadVariableNames',false,'delimiter','tab');
            if size(GG,2) ~= 2
                mess = GenerateMessageEntry(mess, sprintf('ERROR: The globals file ''%s'' contains %g columns, but exactly 2 are expected.',IO.globvar_edit, size(GG,2)));
                flg= true;
            end
            if size(GG,1) ~= nP
                mess = GenerateMessageEntry(mess, sprintf('ERROR: The globals file ''%s'' contains %g rows, but exactly %g are expected to match the number of images provided.',IO.globvar_edit, size(GG,1), nP));
                flg= true;
            end
            PP = table2cell(GG(:,1));
            g = table2array(GG(:,2));
            if ~isnumeric(g)
                mess = GenerateMessageEntry(mess, sprintf('ERROR: The conversion of the 2nd column of the globals file ''%s'' to numeric data failed.',IO.globvar_edit));
                flg= true;
            end
            ind = 1:nP; sPP = cellstr(IO.PP);
            for i=1:nP
                ii = strcmp(sPP,deblank(PP{i,:}));
                iif = find(ii);
                if isempty(iif),
                    mess = GenerateMessageEntry(mess, sprintf('ERROR during matching of global values to images. Could not find case ''%s'' in the IDs of the globals text file!',deblank(IO.PP(i,:))));
                    flg = true;
                    continue
                elseif numel(iif)>1
                    for j=1:numel(iif)
                        mess = GenerateMessageEntry(mess, sprintf('ERROR: Duplicate images detected in position %g: %s!',iif(j),deblank(PP{iif(j)})));
                        flg=true;
                    end
                    continue 
                end
                ind(i)=iif;
            end
         
    end
    if ~flg
        IO.g = g(ind);
        sz_str = sprintf('%s (%g cases): mean (SD) global value: %1.1f (%1.1f)', globvar_str, nP, mean(IO.g), std(IO.g));
    end
    
catch ERR
    mess = GenerateMessageEntry(mess, sprintf('ERROR in %s: \nError message: %s (%s)\nFunction:%s, Line: %g', ...
        'check globals', ERR.message, ERR.identifier, ERR.stack(1).name, ERR.stack(1).line)); 
end

