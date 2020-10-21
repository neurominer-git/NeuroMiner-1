function [Y, Thresh] = nk_ReturnSubSpaces(V, Vm, n_samples, n_subjects, Thresh, reportfile)
% ===============================================================================
% function [Y, Thresh] = nk_ReturnSubSpaces(V, Vm, n_samples, n_subjects, Thresh)
% ===============================================================================
% Vectorizes 3D NIFTI/ANALYZE images based on a space-defining image. This
% is the core image vectorization tool of NM and it currently works only
% together with SPM. 
% 
% (binarized or atlas-based image)
% V             :   Volume structure (output of spm_vol) of images to be read-in
% Vm            :   Volume structure of brain mask / atlas
% n_samples     :   N subject groups
% n_subjects    :   Vector of N x m subjects 
% Thresh        :   Threshold structure with fields: 
%                       Vml (unique values in mask, double), 
%                       nVml (number of unique values in mask, double ), 
%                       Lm (Value descriptors, nVml x 1 cell array), and
%                       threshop (Thresholding operation, cell)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Nikolaoas Koutsouleris, 03/2017

if ~exist('V', 'var') || isempty(V),
    n_samples=1; n_subjects=Inf;
    V = spm_vol(spm_select(n_subjects,'image','Select images to read-in'));
end
if ~exist('Vm', 'var') || isempty(Vm),
    Vm = spm_vol(spm_select(1,'image','Select space-defining image'));
end
if ~exist('Thresh','var') || isempty(Thresh)
    Thresh = nk_DataIO3_SpaceDefImage_config(Vm);
end
if ~exist('reportfile','var'), reportfile = []; end

Lm = Thresh.Lm;
Vml = Thresh.Vml;
nVml = Thresh.nVml;
if nVml~=numel(Vml), nVml=1; end
threshop = Thresh.threshop;

% Prepare data container
Y = cell(1,nVml);

% Check which cell contains image info if V is a cell array
if iscell(V)
    indV=0;
    for i=1:numel(V)
        if ~isempty(V{i})
            indV = i; break
        end
    end
else
    indV=1;
end
if ~indV, error('\nNo images found to read in!'); end

% work in Z direction through the image space
for sl=1:Vm.dim(3)
    
    % Read mask in Z direction by reslicing the subject space to the
    % space-defining image
    M = spm_matrix([0 0 sl 0 0 0 1 1 1]);
    if iscell(V)
        M1  = Vm.mat\V{indV}(1).mat\M;
    else
        M1  = Vm.mat\V(indV).mat\M;
    end
    mask_slice = spm_slice_vol(Vm,M,Vm.dim(1:2),1);
    fprintf('\nReading slice %g/%g',sl,Vm.dim(3))
    
    % Loop though unique labels in the space-defining image
    for l = 1:nVml
        
        % find mask values with current label
        ind0 = find(feval(threshop{1},mask_slice,Vml(l)));        
        % read data inside mask
        if ~isempty(ind0)
            fprintf('\n\tReading label %s (%g voxels): ', Lm{l}, numel(ind0(:)));
            yslice = [];
            
            % Loop through groups of subjects (n_samples)
            for j=1:n_samples
                y = zeros(n_subjects(j), length(ind0));
                fprintf('\n\t\tSample %g: ',j)
                
                % Loop through subjects group j
                for i = 1:n_subjects(j)
                    fprintf('.')
                    try
                        if iscell(V)
                            d = spm_slice_vol( V{j}(i), M1,Vm.dim(1:2), 1 );
                        else
                            d = spm_slice_vol( V(i), M1,Vm.dim(1:2), 1 );
                        end
                    catch
                        error('\nFile problem: %s,',V{j}(i).fname)
                    end
                    y(i,:) = d(ind0);
                end
                yslice = [yslice; y];
            end

            Y{l} = [Y{l} yslice];
        end
    end
    clear mask_slice
end

if ~isempty(reportfile) && exist(reportfile,'file')
    fid = fopen(reportfile,'a');
    fprintf(fid,'\n\nIMPORT SUMMARY:');
    fprintf(fid,'\n===============');
    cLm = char(Lm); sLm = size(cLm,2);
    for l=1:nVml
        fprintf(fid,['\nLabel %' num2str(numel(num2str(nVml))) 'g [ %' num2str(sLm) 's ]: %g voxels'],l, Lm{l}, size(Y{l},2));
    end
end
if nVml <2, Y = Y{1}; end