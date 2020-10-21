function Y = nk_ReadVol2Vec(Vm, V, n_samples, n_subjects, badcoords, verbose)
%
% function Y = nk_ReadVol2Vec(Vm, V, n_samples, n_subjects, verbose)
% 
% This function reads 3D images of a variable number of samples, 
% consisting each of a variable number of subjects into matrix Y
% 
% Inputs:
% -------
% Vm            Space-defining image (brainmask)
% V             Struct Array (SPM) of 3D images to be read in
% n_samples     Number of samples
% n_subjects    Number of subjects in each sample
% badcoords     1xp index to unwanted voxels
% verbose       Verbosity (0 / 1)
% 
% Outputs:
% --------
% Y             m*p matrix with m rows (=observations) and p columns
%               (=voxels), p is eventually reduced due to badcoords.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Nikolaos Koutsouleris, 01/2011

Y=[];
if ~exist('verbose','var') || isempty(verbose)
    verbose = true;
end

if nargin == 0, help nk_ReadVol2Vec, return; end

for sl=1:Vm.dim(3)
    % read mask
    M = spm_matrix([0 0 sl 0 0 0 1 1 1]);
    M1  = Vm.mat\V{1}(1).mat\M;
    mask_slice = spm_slice_vol(Vm,M,Vm.dim(1:2),1);
    if verbose, fprintf('\nReading slice %g/%g: ',sl,Vm.dim(3)); end
    ind0 = find(mask_slice > 0.5);
    clear mask_slice

    % read data inside mask
    if ~isempty(ind0)
        yslice = [];
        for j=1:n_samples
            y = zeros(n_subjects(j), length(ind0));
            if verbose, fprintf(1,'\nsample%g ',j); end
            for i = 1:n_subjects(j)
                if verbose,fprintf('.'), end	
                d = spm_slice_vol(V{j}(i),M1,Vm.dim(1:2),1);
                y(i,:) = d(ind0);
            end
            yslice = [yslice; y];
        end
        Y = [Y yslice];
    end
end

% Remove unwanted voxels
if exist('badcoords','var') && ~isempty(badcoords), Y = Y(:,~badcoords); end

return