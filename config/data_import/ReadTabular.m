function [IO, mess] = ReadTabular(IO, fileflg, mess)

if ~exist('mess','var'), mess=[];end

switch fileflg
    case 1
        IO.Y = evalin('base',IO.M_edit);
        IO.s_featnames = evalin('base',IO.col_edit); 
        IO.s_ID = evalin('base',IO.case_edit);
        if ~isfield(IO,'oocvflag') || ~IO.oocvflag || IO.labels_known
            IO.s_label = evalin('base',IO.label_edit);
        end
    case 2
        load(IO.M_edit);
        IO.Y = eval(IO.matrix_edit);
        IO.s_featnames = eval(IO.col_edit);
        IO.s_ID = eval(IO.case_edit);
        if ~isfield(IO,'oocvflag') || ~IO.oocvflag || IO.labels_known
            IO.s_label = eval(IO.label_edit);
        end
    case {3,4}
        switch fileflg
            case 3
                IO.Y = readtable(IO.M_edit, 'delimiter', IO.delimit);
            case 4
                IO.Y = readtable(IO.M_edit, 'Sheet', IO.sheet);
        end
        IO.s_featnames = IO.Y.Properties.VariableNames;
        col_cases = strcmp(IO.s_featnames,IO.case_edit);
        col_label = RetCellInCellsIndex(IO.s_featnames, IO.label_edit);
        if isnumeric(IO.Y{:,col_cases})
            IO.s_ID = cellstr(num2str(IO.Y{:,col_cases}));
        else
            IO.s_ID = table2cell(IO.Y(:,col_cases));
        end
        if ~isfield(IO,'oocvflag') || ~IO.oocvflag || IO.labels_known
            IO.s_label = table2array(IO.Y(:,col_label));
        end
        IO.Y(:, col_cases | col_label) = []; IO.Y = table2array(IO.Y);
        IO.s_featnames(col_cases | col_label) = [];
end

% Check format of ID variable
if sum(cellfun(@isnumeric,IO.s_ID)) == numel(IO.s_ID)
    IO.s_ID = cellfun(@num2str,IO.s_ID,'UniformOutput',0);
    mess = GenerateMessageEntry(mess, '* Convert IDs: numbers => strings.',0,'blue');
end
if sum(cellfun(@ischar,IO.s_ID)) ~= numel(IO.s_ID)
    mess = GenerateMessageEntry(mess, 'ERROR during data import: The ID variable should contain only alpha(numeric) strings or numbers.');
    return 
end

% Check whether previous feature names exists and make sure that the order
% of features in the new modality matches the order of features in previous modalities
if isfield(IO,'oocvflag') && IO.oocvflag && isfield(IO,'featnames')
    if ~IO.feat_reordered
        ind = zeros(1,numel(IO.featnames)); flg = false;
        for i=1:numel(IO.featnames)
            ii = strcmp(IO.s_featnames,IO.featnames{i});
            if ~any(ii)
                mess = GenerateMessageEntry(mess, 'Error during import of new matrix-based modality data.\n\tCould not find feature ''%s'' among the features of the new modality!');
                flg = true;
                continue;
            end
            ind(i) = find(ii);
        end
        if flg, return; end
        if ~isequal(ind,1:numel(IO.featnames))
            IO.Y = IO.Y(:,ind);
            IO.features_reordered = ind;
            mess = GenerateMessageEntry(mess, sprintf('Data (%s) reordered during import to match previous feature order!',IO.desc_data), 0, 'blue'); 
        end
    end
else
    IO.featnames = IO.s_featnames;
end

% Check whether previous label exists and make sure that the order of cases
% in the new modality matches the order of cases in previous modalities

if isfield(IO,'ID') 
    ind = zeros(numel(IO.ID),1); flg = false; 
    for i=1:numel(IO.ID)
        ii = find(strcmp(IO.s_ID,IO.ID{i}));
        if ~any(ii),
            mess = GenerateMessageEntry(mess, sprintf('Error during import of new matrix-based modality data. Could not find case ''%s'' in the IDs of the new modality!',IO.ID{i}));
            flg = true;
            continue
        elseif numel(ii)>1
            mess = GenerateMessageEntry(mess, sprintf('Error during import of new matrix-based modality data. Found %g instances of case ''%s'' in the IDs of the new modality!', numel(ii), IO.ID{i}));
            flg = true; 
            continue
        end
        ind(i) = ii; 
    end
    if flg, return; end
    if ~isequal(ind,1:numel(IO.ID))
        IO.Y = IO.Y(uint16(ind),:);
        IO.cases_reordered=ind;
        IO.ID = IO.s_ID(ind);
        mess = GenerateMessageEntry(mess, sprintf('Data (%s) reordered during import to match previous case order!',IO.desc_data),0,'blue'); 
    end

else
    IO.ID = IO.s_ID;
    if ~isfield(IO,'oocvflag') || ~IO.oocvflag || IO.labels_known
        switch IO.modeflag
            case 'classification'           
                nL = size(IO.s_label,2);
                IO.label = nan(size(IO.s_label));     
                if iscell(IO.s_label)
                    IO.nonnan_label = ~(isnanincell(IO.s_label) | strcmp(IO.s_label,'') | cellfun(@isempty,IO.s_label));
                else
                    IO.nonnan_label = ~isnan(IO.s_label);
                end
                if nL>1
                    IO.n_subjects = cell(nL,1);
                else
                    IO.n_subjects = zeros(1,numel(unique(IO.s_label(IO.nonnan_label),'stable'))); 
                end
                tdesc_groups = IO.desc_groups;
                for j=1:nL
                    desc_groups = unique(IO.s_label(IO.nonnan_label(:,j),j),'stable');
                    IO.n_samples(j) = numel(desc_groups);
                    if nL>1
                        IO.desc_groups{j} = desc_groups;
                    else
                        IO.desc_groups = desc_groups;
                    end
                    for i=1:numel(tdesc_groups)
                        if isfield(IO,'selCases') && ~isempty(IO.selCases)
                            selCases = IO.selCases;
                        else
                            selCases = true(size(IO.s_label,1),1);
                        end
                        ind = strcmp(IO.s_label(:,j),tdesc_groups{i}) & selCases;
                        IO.label(ind,j) = i;
                        if nL>1
                            IO.n_subjects{j}(i) = sum(ind); 
                        else
                            IO.n_subjects(i) = sum(ind); 
                        end
                    end
                end
            case 'regression'
                nL = size(IO.s_label,2);
                for j=1:nL
                    if ~isfield(IO,'desc_groups') || numel(IO.desc_groups)<nL
                        IO.desc_groups{j} = nk_input('Describe the sample used to construct the regression model',0,'s');
                    end
                end
                if isfield(IO,'selCases') && ~isempty(IO.selCases)
                    selCases = IO.selCases;
                else
                    selCases = true(size(IO.s_label,1),1);
                end
                IO.label = IO.s_label(selCases);
                IO.n_subjects = size(IO.label,1);
        end
    end
end

IO.n_subjects_all = sum(IO.n_subjects);