function [check_str, mess, t_label, IO] = CheckLabel(IO, act, t_Y, na_str, mess, Mx)

check_str = na_str;
t_label = [];
if ~isfield(IO,'label_edit') || any(strcmp(IO.label_edit,na_str)), return, end
if ~exist('mess','var'), mess=[]; end
if ~exist('Mx','var') || isempty(Mx), Mx = 'unknown matrix source'; end
if ~exist('t_Y','var') || isempty(t_Y), 
    mess = GenerateMessageEntry(mess, sprintf('ERROR in %s: Input data is missing!', act)); 
    return
else
    t_Y_sz = size(t_Y);
end
if ~iscell(IO.label_edit)
    label_edit = cellstr(IO.label_edit);
    nL = 1;
else
    label_edit = IO.label_edit;
    nL = numel(IO.label_edit);
end

try
    for j=1:nL
        jn_subjects = []; jdesc_groups = [];
        if isfield(IO,'n_subjects')
            if iscell(IO.n_subjects)
                jn_subjects = IO.n_subjects{j};
                jdesc_groups = IO.desc_groups{j};
            else
                jn_subjects = IO.n_subjects;
                jdesc_groups = IO.desc_groups;
            end
        end
        switch IO.groupmode
            case {0, 1}
                t_label = evalin('base', label_edit{j}); 
            case 2
                load(IO.M_edit, label_edit{j}); t_label = eval(label_edit{j}); 
            case {3,4}
                t_label = table2array(t_Y(:,label_edit(j)));
        end
        
        % Ceck whether unlabeled datapoints have to be added to the label
        % matrix
        if IO.nangroup == 1 && IO.nan_subjects > 0
           t_label = [t_label; nan(IO.nan_subjects,size(t_label,2))];
        end
        
        % Check whether Nan is included in the data
        if iscell(t_label)
            ind_nan = isnanincell(t_label) | strcmp(t_label,'');
        else 
            ind_nan = logical(sum(isnan(t_label),2));
        end
        IO.nonnan_label = ~ind_nan;

        switch IO.modeflag
            case 'classification'
                if ~iscell(t_label), 
                    mess = GenerateMessageEntry(mess, sprintf('ERROR in %s: Var ''%s'' is not a cell array of strings!', act, IO.label_edit));
                    return; 
                end    
            case 'regression'
                if ~isnumeric(t_label), 
                    mess = GenerateMessageEntry(mess, sprintf('ERROR in %s: Var ''%s'' is not a numeric array!', act, IO.label_edit)); 
                    return; 
                end
        end

        if ~isempty(t_Y) && size(t_label,1) ~= t_Y_sz(1)
            mess = GenerateMessageEntry(mess, sprintf('ERROR in %s: Number of cases in label array ''%s'' does not match number of cases (%g) in the data!', act, label_edit{j}, t_Y_sz(1))); 
            return;  
        end
        
        if isfield(IO,'selCases') && ~isempty(IO.selCases)
            selCases = IO.selCases;
        else
            selCases = true(size(t_label,1),1);
        end

        if strcmp(IO.modeflag,'classification') 
            desc_groups = unique(t_label(IO.nonnan_label),'stable');
            IO.desc_groups = desc_groups;
            n_samples = numel(desc_groups);
            for i=1:n_samples
                n_subjects(i) = sum(strcmp(t_label(IO.nonnan_label),desc_groups{i}) & selCases);
            end
        else
            desc_groups = IO.desc_groups;
            n_samples = 1; n_subjects = size(t_label(IO.nonnan_label & selCases,:),1);
        end
        
        % Perform a series of label checks if previous modalities exist
        if ~isempty(jn_subjects) && ~any(~isfinite(jn_subjects))
           % Check whether the number of encoded samples matches
           if n_samples ~= numel(jn_subjects)
             mess = GenerateMessageEntry(mess, sprintf('ERROR in %s: Number of samples (%g) encoded in ''%s'' does not match the number of samples in previous modalities (%g)!', ...
                act, n_samples, label_edit{j}, numel(jn_subjects))); 
             return;  
           end
           flg = false;
           % Check whether label names and sample sizes match by looping through each sample
           for i=1:n_samples
               %jdesc_groups_i = find(strcmp(jdesc_groups, desc_groups{i}));
               if ~any(strcmp(jdesc_groups, desc_groups{i}))
                    mess = GenerateMessageEntry(mess, sprintf('ERROR in %s: Name (%s) of group %g does not match any group name in pre-existing label!', ...
                        act, desc_groups{i}, i)); 
                    flg = true;
               end
%                if n_subjects(i) ~= jn_subjects(i) 
%                    mess = GenerateMessageEntry(mess, sprintf('ERROR in %s: Number of cases (%g) assigned to label ''%s'' does not match the pre-existing number of cases for this label (%g)!', ...
%                         act, n_subjects(i), jdesc_groups{i}, jn_subjects(i))); 
%                     flg = true;
%                end
               % Check whether the cases of group i in previous modality
               % match cases of group i in current modality
               if strcmp(IO.modeflag,'classification') && ~strcmp(IO.case_edit,na_str)
                   II = strcmp(t_label,desc_groups{i});
                   switch IO.groupmode
                       case {0,1}
                           curr_iID = evalin('base', IO.case_edit); 
                           curr_iID = curr_iID(II); 
                       case 2
                           curr_iID = load(IO.M_edit, IO.case_edit); 
                           curr_iID = curr_iID(II); 
                       case {3,4}
                           IDcol = strcmp(t_Y.Properties.VariableNames, IO.case_edit);
                           if isnumeric(t_Y{II,IDcol})
                               curr_iID = cellstr(num2str(t_Y{II,IDcol}));
                           else
                               curr_iID = cellstr(t_Y{II,IDcol}) ;
                           end
                   end
                   if numel(IO.ID) == size(t_label,1),
                       prev_iID = IO.ID(IO.label==find(strcmp(jdesc_groups,desc_groups{i})));
                       MIS = setdiff(curr_iID, prev_iID);
                       if ~isempty(MIS)
                           mis_str = sprintf('\n=> %s',strjoin(MIS,'\n=> '));
                           mess = GenerateMessageEntry(mess, sprintf('ERROR in %s: Cases of sample %g (''%s'') do not match cases previously assigned to this sample:%s', act, i, jdesc_groups{i},mis_str)); 
                           flg = true;
                       end
                   end
               end
           end
           if flg, return; end
        end
        
        if any(ind_nan)
            mess = GenerateMessageEntry(mess, sprintf(['WARNING in %s: label array ''%s'' contains %g label(s) with non-finite entries!\n' ...
                                     'Consider activating label imputation methods in the NM workspace setup'], act, IO.label_edit, sum(ind_nan)), 2, 'SystemCommands'); 
        end

        
    end

catch ERR
    t_label=[];
    mess = GenerateMessageEntry(mess, sprintf('ERROR in %s: Label data ''%s'' could not be extracted from %s.\nError message: %s (%s)\nFunction:%s, Line: %g', ...
        act, label_edit{j}, Mx, ERR.message, ERR.identifier, ERR.stack(1).name, ERR.stack(1).line)); 
    return
end
check_str = '';
if strcmp(IO.modeflag,'classification')
    for j=1:nL
        if nL>1
            check_str = sprintf('%s; (%g) Var ''%s'': %g groups', check_str, j, label_edit{j}, n_samples);
        else
            check_str = sprintf('%s; Var ''%s'': %g groups', check_str, label_edit{j}, n_samples);
        end
        for i=1:n_samples
            check_str = sprintf('%s, %s (N=%g)', check_str, desc_groups{i}, n_subjects(i));
        end
    end
else
    for j=1:nL
        if size(t_label,2)>1
            check_str = sprintf('%s; Var ''%s'': %g cases with %g labels (%g unique values)', check_str, label_edit{j}, n_subjects, size(t_label,2), numel(unique(t_label(isfinite(t_label)))));
                   else
            check_str = sprintf('%s; Var ''%s'': %g cases with %g unique values', check_str, label_edit{j}, n_subjects, numel(unique(t_label(isfinite(t_label)))));
        end
    end
end

if any(~IO.nonnan_label), 
    check_str = sprintf('%s; %g undefined label values(s)',check_str, sum(~IO.nonnan_label)); 
end

check_str(1:2)=[];

mess = GenerateMessageEntry(mess, '* Label tests passed', 0, 'blue');