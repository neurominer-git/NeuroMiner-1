function [PP, P, V, t_subjects, n_dims, vox] = nk_DefineFiles(dtype, n_subjects, n_samples, groupnames, filt, labelflag, addflag, n_dims, vox)

PP=[]; P = cell(n_samples); V = cell(n_samples, 1);

if ~exist('n_dims','var'), n_dims = []; end
if ~exist('vox','var'), vox = []; end
t_subjects = zeros(n_samples,1);

% Loop through samples and invoke file selection dialogue
for i = 1:n_samples
    if labelflag
        hdrstr = ['Select images for group ' num2str(i) ': ' groupnames{i}];
        if addflag
            nsubj = Inf;
        else
            nsubj = n_subjects(i);
        end
    else
        hdrstr = 'Select images';
        nsubj = Inf;
    end
    [P{i}, V{i}] = nk_FileSelector(nsubj, dtype, hdrstr, filt);
    if isempty(P{i}), 
        error('No images specified for group %s!', groupnames{i}); 
    end
    t_subjects(i) = size(P{i},1);
    PP = char(PP,P{i});	
end

[errfl, n_dims, vox] = check_image_dimensions(V, n_dims, vox);
switch errfl
    case 0
        fprintf('\n');cprintf('green','Image checks passed! All images have the same size and voxel size')
    case 1
        error('Incompatible image dimensions')
    case 2
        error('Incompatible voxel sizes')
end
PP = PP(2:end,:);

function [errfl, n_dims, vox] = check_image_dimensions(V, n_dims, vox)

n_samples = numel(V); errfl = 0;
for i = 1:n_samples
    for j=1:numel(V{i});
        
        % Check image dimensions
        if isempty(n_dims), 
            n_dims = prod(V{i}(j).dim); 
        elseif n_dims ~= prod(V{i}(j).dim)
            errfl = 1; break
        end
        
        % Check voxel sizes
        if isempty(vox),
            vox = diag(V{i}(j).mat);
        elseif ~isequal(vox, diag(V{i}(j).mat))
             errfl = 2; break
        end
    end
end