function [sz_str, mess, IO] = CheckTabFile(IO, act, na_str, mess)

if ~exist('mess','var'), mess=[]; end
sz_str = na_str;

warning('off')
try
    switch IO.groupmode
        case 1
            Mx = 'MATLAB workspace'; matrixname = IO.M_edit;
           
        case 2
            Mx = 'MATLAB file'; matrixname = IO.matrix_edit;
           
        case 3
            Mx = 'text file'; matrixname = IO.M_edit;
           
        case 4
            Mx = 'spreadsheet file'; matrixname = IO.M_edit;
    end
    if (strcmp(act,'check_matrix') && (IO.groupmode == 3 || IO.groupmode == 4)) && (isfield(IO,'oocvflag') && IO.oocvflag )
       unknown_str=''; if ~IO.labels_known, unknown_str = ' (unused because you specified that the labels of the new data are unknown)'; end
       mess = GenerateMessageEntry(mess, sprintf(['INFO: The ''%s'' containing the independent test data should be defined as follows:' ...
                                 '\n\t* Case IDs in a column named ''%s''', ...
                                 '\n\t* Labels in a column named ''%s''%s', ...
                                 '\n\t* %g feature columns with the same names as in the discovery %s'], ...
                                 matrixname, IO.case_edit, IO.label_edit, unknown_str, numel(IO.featnames), Mx),0,'blue');
    end
    source = Mx;
catch ERR
    mess = GenerateMessageEntry(mess, sprintf('ERROR in %s: \nError message: %s (%s)\nFunction:%s, Line: %g', ...
        act, ERR.message, ERR.identifier, ERR.stack(1).name, ERR.stack(1).line)); 
    return
end

if ~isfield(IO,'t_Y') || isempty(IO.t_Y) || IO.reread_mat
    try
        switch IO.groupmode
            case 1
                if ~strcmp(IO.M_edit,na_str), 
                    IO.t_Y = evalin('base', IO.M_edit);
                else
                    return
                end
                if ~strcmp(IO.col_edit,na_str), 
                    try
                        s_featnames = evalin('base',IO.col_edit); 
                    catch
                        mess = GenerateMessageEntry(mess, sprintf('ERROR in %s: Feature name vector ''%s'' not found in %s!', act, IO.col_edit, source));
                        return
                    end
                end
            case 2
                load(IO.M_edit, IO.matrix_edit); 
                IO.t_Y = eval(IO.matrix_edit);
                try
                    if ~strcmp(IO.col_edit,na_str), s_featnames = eval(IO.col_edit); end
                catch
                    mess = GenerateMessageEntry(mess, sprintf('ERROR in %s: Feature name vector ''%s'' not found in %s!', act, IO.col_edit, source));
                    return
                end
            case 3
                IO.t_Y = readtable(IO.M_edit, 'Delimiter', IO.delimit);
                s_featnames = IO.t_Y.Properties.VariableNames;
            case 4
                if strcmp(IO.sheet,na_str)
                    mess = GenerateMessageEntry(mess, sprintf('ERROR %s: You have not specified a source sheet in ''%s''.',act, matrixname)); 
                    return
                else
                    IO.t_Y = readtable(IO.M_edit, 'Sheet', char(IO.sheet));
                    s_featnames = IO.t_Y.Properties.VariableNames;
                    fprintf('\nChecking spreadsheet file %s ...',IO.M_edit);
                    [~,IO.sheets] = xlsfinfo(IO.M_edit);
                end
        end
    catch
        switch IO.groupmode
            case {1,2}
                mess = GenerateMessageEntry(mess, sprintf('ERROR in %s: Matrix ''%s'' not found in %s!', act, matrixname, source));
            otherwise
                mess = GenerateMessageEntry(mess, sprintf('ERROR in %s: Could not open ''%s''.\n\tPlease check whether the file exists, is corrupt or the import settings are not completed.', act,IO.M_edit)); 
        end
        IO.reread_mat = true;
        return
    end
    if IO.oocvflag 
        if IO.labels_known, Li = 2; else, Li=[]; end
        ind_feat=nan(1,numel(IO.featnames_cv)); notfnd = false;
        for f=1:numel(IO.featnames_cv)
            fI = find(strcmp(s_featnames,IO.featnames_cv{f}));
            if ~isempty(fI), 
                ind_feat(f) = fI;
            else
                mess = GenerateMessageEntry(mess, sprintf('Could not find feature ''%s'' of the discovery data in new data matrix.',IO.featnames_cv{f})); 
                notfnd = true;
            end
        end
        if notfnd 
            IO.reread_mat = true;
            return
        else
            switch IO.groupmode
                case {1,2}
                    IO.t_Y = IO.t_Y(:,ind_feat);
                case {3,4}
                    IO.t_Y = [IO.t_Y(:,1) IO.t_Y(:,Li) IO.t_Y(:,ind_feat)];
            end
           
            IO.featnames = s_featnames(ind_feat);
            IO.selFeats = true(size(IO.featnames));
            IO.feat_reordered = true;
            mess = GenerateMessageEntry(mess, 'Successfully extracted features corresponding to discovery matrix!',0,'blue'); 
        end
        
    end
    if size(IO.t_Y,2) < 2
        switch IO.groupmode
            case {1,2}
                mess = GenerateMessageEntry(mess, sprintf('WARNING in %s: I found only one column in the ''%s''.\n\tPlease check whether you specified the correct matrix.', act, source), 2, 'SystemCommands');
            case {3,4}
                mess = GenerateMessageEntry(mess, sprintf('WARNING in %s: I found only one column in the ''%s''.\n\tPlease check whether the delimiter (%s) is correct.', act, source, IO.delimit), 2,'SystemCommands');
        end
        IO.reread_mat = true;
    else
        IO.reread_mat = false;
    end
end

t_Y_sz = size(IO.t_Y);

try
    switch act  
    
    case 'check_label'
        [sz_str, mess, ~, IO] = CheckLabel(IO, act, IO.t_Y, na_str, mess, Mx);
        
    case 'check_ID'
        if ~isfield(IO,'case_edit') || strcmp(IO.case_edit,na_str), return, end
        try
            switch IO.groupmode
                case 1
                    t_ID = evalin('base',IO.case_edit); 
                case 2
                    load(IO.M_edit, IO.case_edit); t_ID =eval(IO.case_edit);
                case {3,4}
                    t_ID = table2cell(IO.t_Y(:,{IO.case_edit}));
            end
            
            t_ID_sz = size(t_ID);
            if ~ismatrix(t_ID_sz) || size(t_ID,2)>1 , 
                mess = GenerateMessageEntry(mess, sprintf('ERROR in %s: cases vector should not have more than one dimension!', act));
                return
            end
            if ~iscell(t_ID) || sum(cellfun(@ischar,t_ID))~=numel(t_ID)
                mess = GenerateMessageEntry(mess, sprintf('Warning in %s: cases vector is not a cell array of strings or has non-string elements! Conversion will be performed during import.', act),2,'SystemCommands');
            end
            if ~isempty(IO.t_Y) && t_ID_sz(1) ~= t_Y_sz(1)
                mess = GenerateMessageEntry(mess, sprintf('ERROR in %s: Number of cases in cases array ''%s'' does not match number of cases (%g) in ''%s''!', act, IO.case_edit, t_Y_sz(1), matrixname));
                return
            end
            sz_str = sprintf('Var ''%s'': %g cases',IO.case_edit, numel(t_ID));
        catch
            mess = GenerateMessageEntry(mess, sprintf('ERROR in %s: Case ID variable(s) ''%s'' not found in %s.', act, IO.case_edit, Mx)); 
            return
        end
        mess = GenerateMessageEntry(mess, '* ID tests passed', 0, 'blue');
        
    case 'check_featnames'
        if ~isfield(IO,'col_edit') || strcmp(IO.col_edit,na_str), return, end
        try
            switch IO.groupmode
                case 1
                    t_featnames = evalin('base',IO.col_edit); 
                case 2
                    load(IO.M_edit, IO.col_edit); t_featnames = eval(IO.col_edit); 
            end
            t_featnames_sz = size(t_featnames);
            if numel(t_featnames_sz)>2 || ~any(t_featnames_sz == 1) , 
                mess = GenerateMessageEntry(mess, sprintf('ERROR %s: feature names vector should not have more than one dimension!', act)); 
                return
            end
            if ~iscell(t_featnames), 
                mess = GenerateMessageEntry(mess, sprintf('ERROR %s: feature names vector is not a cell array of strings!', act)); 
                return
            end
            if ~isempty(t_featnames) && numel(t_featnames) ~= t_Y_sz(2)
                mess = GenerateMessageEntry(mess, sprintf('ERROR in %s: Number of features in feature names array ''%s'' does not match number of features (%g) in ''%s''!', act, IO.col_edit, t_Y_sz(2), source)); 
                return
            end
            sz_str = sprintf('Var ''%s'': %g features',IO.col_edit,numel(t_featnames));
        catch
            mess = GenerateMessageEntry(mess, sprintf('Error in %s: Feature names variable ''%s'' not found in ''%s''!', act, IO.col_edit, source));
            return
        end
        mess = GenerateMessageEntry(mess, '* Feature name tests passed', 0, 'blue');
        
    case 'check_matrix'
        
        switch IO.groupmode
            case 1
                if ~isfield(IO,'M_edit'), return, end
                IO.t_Y = evalin('base',IO.M_edit);
                t_Y = IO.t_Y;
            case 2
                if ~isfield(IO,'matrix_edit'), return, end
                load(IO.M_edit,IO.matrix_edit); IO.t_Y =eval(IO.matrix_edit); 
                t_Y = IO.t_Y;
            case {3,4}
                if ~isfield(IO,'M_edit'), return, end
                if ~isfield(IO,'case_edit'), return, end
                if ~isfield(IO,'label_edit'), return, end
                featnames = IO.t_Y.Properties.VariableNames;
                col_cases = strcmp(featnames,IO.case_edit);
                col_label = RetCellInCellsIndex(featnames, IO.label_edit);
                t_Y = table2array(IO.t_Y(:, ~col_cases & ~col_label ));
        end

        if ~ismatrix(t_Y_sz),
            mess = GenerateMessageEntry(mess, sprintf('ERROR in %s: Input matrix should not have more than two dimensions!', act));
            return
        end

        if iscell(t_Y)
            mess = GenerateMessageEntry(mess, sprintf('ERROR  in %s: Matrix is a cell array! Check for non-numeric values in the source data.', act)); 
            return
        end
        
        sumNonNum = sum(~isfinite(t_Y(:)));
        if sumNonNum 
            percNonNum = sumNonNum*100/numel(t_Y(:));
            mess = GenerateMessageEntry(mess, ...
                sprintf('WARNING in %s: Matrix contains %g (%1.1f percent of data) non-finite value(s).\n\tUse imputation methods if data are OK.', act, sumNonNum, percNonNum), 2, 'SystemCommands');
        end

        sz_str = sprintf('Matrix in ''%s'': %g cases, %g features', matrixname, t_Y_sz(1),t_Y_sz(2));
        
        if isfield(IO,'selFeats') && ~isempty(IO.selFeats) && sum(IO.selFeats)<numel(IO.selFeats)
            mess = GenerateMessageEntry(mess, sprintf('* %g (of %g) features will be removed from the matrix based on your selection.', ...
                numel(IO.selFeats) -sum(IO.selFeats), numel(IO.selFeats) ),0,'blue'); 
        end
        if isfield(IO,'selCases') && ~isempty(IO.selCases) && sum(IO.selCases)<numel(IO.selCases)
            mess = GenerateMessageEntry(mess, sprintf('* %g (of %g) cases will be removed from the matrix based on your selection.', ...
                numel(IO.selCases) -sum(IO.selCases), numel(IO.selCases) ),0,'blue'); 
        end
        
        mess = GenerateMessageEntry(mess, '* Predictor matrix tests passed', 0, 'blue');
    end
    
    
    
catch ERR
    mess = GenerateMessageEntry(mess, sprintf('\n ERROR caused by action %s:\n Message: %s (%s)\n Function: %s\n Line: %g', act, ERR.message, ERR.identifier, ERR.stack(1).name, ERR.stack(1).line));
    return
end

