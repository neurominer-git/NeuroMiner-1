function [n_samples, n_subjects, groupnames, labelflag] = nk_DefineGroups(oocv, modeflag, n_subjects, groupnames)

switch modeflag
    case 'classification'
        if (~exist('groupnames','var') || isempty(groupnames)) 
            n_samples = nk_input('# of groups?',0,'e');
        else
            n_samples = numel(groupnames);
        end
        if (~exist('n_subjects','var') || isempty(n_subjects)) 
            n_subjects = Inf(1,n_samples);
        end
    case 'regression'
        n_samples = 1;
        if (~exist('n_subjects','var') || isempty(n_subjects))
            n_subjects = Inf; 
        end
end
   
if ~exist('groupnames','var') || isempty(groupnames)
    for i=1:n_samples
        if n_samples > 1
            groupnames{i} = nk_input(sprintf('Define name of group #%g',i),0,'s');
        else
            groupnames{i} = nk_input( 'Define name of sample',0,'s');
        end
    end
end

switch oocv
    case 0
       labelflag = true;
    case 1
       labelflag = nk_input('Do you know the labels of these cases',0,'yes|no',[1,0],1);
       if strcmp(modeflag,'classification')   
           if ~labelflag, n_samples = 1; end
       else
           n_samples = 1;
       end
end