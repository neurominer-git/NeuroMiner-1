function IO = DefineLabels(IO, modeflag)

switch modeflag
    case 'classification'
        if isfield(IO,'L') && ~isempty(IO)
            IO.label = zeros(size(IO.L));
            IO.n_label = size(IO.L,2); 
            uL = unique(IO.L,'stable');
            nLall = numel(uL); 
            for j=1:IO.n_label
                for i=1:nLall
                    if isnumeric(IO.L)  
                        ind = IO.L == uL(i);
                        IO.desc_groups{i} = nk_input(sprintf('Define name of group #%g (N=%g):',i ,sum(ind)),'s');
                    else
                        ind = strcmp(IO.L,uL{i});
                        IO.desc_groups{i} = uL{i};
                    end
                    IO.label(ind,j) = i; 
                end
            end
        elseif isfield(IO,'n_subjects') && isfield(IO,'n_samples') 
            IO.label = [];
            for i=1:IO.n_samples
                IO.label = [ IO.label; i*ones(IO.n_subjects(i),1) ];
            end
        end
    case 'regression'
        if isfield(IO,'L') && ~isfield(IO,'label')
            IO.label = IO.L;
        end
end
% add missing labels as NaNs to the label vector/matrix
if isfield(IO,'nan_subjects') && IO.nan_subjects>0
   IO.label = [IO.label; nan(IO.nan_subjects,size(IO.label,2))]; 
end