function [cases, files, mess] = nk_DefineCaseNames2(P, N, cases_comp, mess) 

if ~exist('mess','var') || isempty(mess), mess=[];end

if ~exist('P','var') || isempty(P)
    cprintf('red','\nOperating in matrix mode.')
    caseopt= nk_input('Where do I find the Subject IDs ?',0,'m', ...
        'In the MATLAB(R) workspace|In a text file|Autogenerate subject IDs',1:3,1);	

    switch caseopt

        case 1
            P = nk_input('Provide MATLAB variable with subject IDs',0,'e',[],sampsize);
            if isnumeric(P), P = num2str(P,'%g'); end
            casemanip = nk_input('extract variable cases substrings',0,'yes|no',[1 0],1);

        case 2
            P = nk_FileSelector(1,'matrix','Select text file with subject IDs','.*\.txt$|.*\.csv');
            P = textscan(P,'%c');
            if ~exist('N','var') || isempty(N), N = size(P,1); end
            if size(P,1) ~= N
                errordlg('Length of subject ID vector does not match length of data matrix!')
            end
            casemanip = nk_input('Shall I try to extract unique substrings from subject IDs',0,'yes|no',[1 0],1);
        case 3
            P = cell(N,1);
            for i = 1:N
                P{i} = sprintf('Case_%1.0f',i);
            end
            casemanip = 0;
    end
    files = P;
else
    cprintf('red','\nOperating in file mode.')
    casemanip = 1;
    files = P;
end

% Perform CaseID extraction
if casemanip
    P = cellstr(P);
    x = spm_file(P,'basename'); nx = numel(unique(x)); nP = size(P,1);
    if nx == nP
        % The filenames contain correct number of unique IDs
        if nP>1
            cases = x; [~,C] = spm_str_manip(cases,'C'); cases = C.m;
        else
            cases = x;
        end
    elseif nx > 1 && nx < nP
        % Duplicates?
        num_occur = arrayfun( @(i) sum(strcmp(x, x{i})), 1:nP);
        ind_num_occur = find(num_occur>1);
        cases_notunique = x(num_occur>1); cases_notunique_str=sprintf('\n\nFollowing cases were repeatedly found in the list:');
        for i=1:numel(cases_notunique)
            cases_notunique_str = sprintf('%s\n\t%s: Position in list: %3g | Number of occurences: %3g', cases_notunique_str,cases_notunique{i}, ind_num_occur(i), num_occur(i));
        end
        mess = GenerateMessageEntry(mess,sprintf(['The number of unique case IDs (N=%g) is smaller than the number of files (N=%g)!\n' ...
                                                  'If case IDs reside in the filenames check whether you entered duplicate filenames.\n' ...
                                                  'If case IDs reside in the paths, but not in the filenames, make sure that this is really the case (filenames should be identical!):%s'],nx,nP,cases_notunique_str));
        cases=[];
        return
    elseif nx == 1
        notfnd = true;
        y = spm_file(P,'path');
        while notfnd  
            x = spm_file(y,'basename');
            if numel(unique(x)) == size(P,1)
                cases = x; [~,C] = spm_str_manip(cases,'C'); cases = C.m; notfnd = false;
            elseif sum(cellfun(@isempty,x)) == size(P,1)
                cases = [];
                mess = GenerateMessageEntry(mess,'Case ID definition: No unique case IDs could be extracted from file paths');
                return
            else
                y = spm_file(y,'path');
            end
        end
    else
        cases = [];
        mess = GenerateMessageEntry(mess,'Case ID definition: No unique case IDs could be extracted from file paths');
        return
    end
   
else
    cases = P;
end

% Check against existing IDs
if exist('cases_comp','var') && ~isempty(cases_comp)
    
    ind = zeros(numel(cases_comp),1);
   
    for i=1:numel(cases)
        ii = strcmp(cases,cases_comp{i});
        if ~any(ii),
            mess = GenerateMessageEntry(mess, sprintf('Warning: Could not find case ''%s'' in the IDs of the new modality!',cases_comp{i}),'SystemCommands', 2); 
            cnt=cnt+1;
            continue
        end
        try
            ind(i) = find(ii);
        catch
            mess = GenerateMessageEntry(mess, sprintf('ERROR: found %g instances of case ''%s'' in the IDs!',numel(find(ii)),cases_comp{i})); 
            return
        end    
    end
    
end
