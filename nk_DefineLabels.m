function [label, n_subjects]= nk_DefineLabels(V, label, n_subjects, modeflag, addfl, CVstr, NaNflag)

fprintf('\n\n*** Label definition ***\n')

if exist('label','var') && ~isempty(label)
    fprintf('\nPrevious %s-subjects'' label vector detected.', CVstr)
    if addfl
       fprintf('\nAdding new labels.')
    else
       return
    end
end

if ~exist('NaNflag','var') && ~isempty(NaNflag), NaNflag = false; end
    
switch modeflag 
    
    case 'classification'
        
        % Automatically create target variables for classification
        fprintf('Automatic label definition [classification].\n')
        tlabel = [];
        for j=1:numel(V)
            n_subjects(j) = size(V{j},1);
            if j==numel(V) && NaNflag
                tlabel = [tlabel; nan(n_subjects(j),1)];
            else
                tlabel = [tlabel; j*ones(n_subjects(j),1)];
            end
        end
        if exist('label','var') && ~isempty(label) && addfl
            label = [label; tlabel];
        else
            label = tlabel;
        end
        
    case 'regression';
        
        if exist('label','var') && ~isempty(label) && addfl
            [n_label, m_label] = size(label);
            Prmpt = ['Add new target labels to existing ones [regression (' CVstr ')]'];
        else
            n_label = 0; label = []; m_label = 1;
            Prmpt = ['Define target labels [regression (' CVstr ')]'];
        end
        
        % Read-in target variables for regression 
        tlabel = nk_input(Prmpt,0,'r',[],[n_subjects-n_label m_label]);
        label = [label; tlabel];
        
        if NaNflag, label = [label; nan(n_subjects(end),1)]; end
end