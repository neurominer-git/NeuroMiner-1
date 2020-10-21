function [T, VO] = nk_ResizeImg(P, Voxdim, BB, FileName)
% nk_ResizeImg -- resample images to have specified voxel dims and BBox
% nk_ResizeImg(P, voxdim, bb)
%
% Output images will be prefixed with 'r', and will have voxel dimensions
% equal to voxdim. Use NaNs to determine voxdims from transformation matrix
% of input image(s).
% If bb == nan(2,3), bounding box will include entire original image
% Origin will move appropriately.

% Based on John Ashburner's reorient.m
% http://www.sph.umich.edu/~nichols/JohnsGems.html#Gem7
% http://www.sph.umich.edu/~nichols/JohnsGems5.html#Gem2

% reslice images one-by-one
if ischar(P), V = spm_vol(P); else V=P; end

% default voxdim to current volume's voxdim, (from mat parameters)
if any(isnan(Voxdim))
    vprm = spm_imatrix(V.mat);
    vvoxdim = vprm(7:9);
    Voxdim(isnan(Voxdim)) = vvoxdim(isnan(Voxdim));
end
Voxdim = abs(Voxdim(:))'; % (handle analyze_flip separately)

bbmn = BB(1,:);
bbmx = BB(2,:);
% default BB to current volume's
if any(isnan(BB(:)))
    d = V.dim(1:3);
    % corners in voxel-space
    c = [ 1    1    1    1
          1    1    d(3) 1
          1    d(2) 1    1
          1    d(2) d(3) 1
          d(1) 1    1    1
          d(1) 1    d(3) 1
          d(1) d(2) 1    1
          d(1) d(2) d(3) 1 ]';
    % corners in world-space
    tc = V.mat(1:3,1:4)*c;
    % reflect in x if required
    if spm_flip_analyze_images; tc(1,:) = -tc(1,:); end;
    % bounding box (world) min and max
    mn = min(tc,[],2)'; % (don't think it's necessary to round these...)
    mx = max(tc,[],2)';
    bbmn(isnan(bbmn)) = mn(isnan(bbmn));
    bbmx(isnan(bbmx)) = mx(isnan(bbmx));
end

% voxel [1 1 1] of output should map to BB mn
% (the combination of matrices below first maps [1 1 1] to [0 0 0])
mat = spm_matrix([bbmn 0 0 0 Voxdim])*spm_matrix([-1 -1 -1]);
% voxel-coords of BB mx gives number of voxels required (now round up)
imgdim = ceil(mat \ [bbmx 1]')';

% reflect in x if required
if spm_flip_analyze_images; mat = diag([-1 1 1 1])*mat; end;

% Create output image
VO            = V;
VO.fname      = fullfile(pwd,'temp.nii');
VO.dim(1:3)   = imgdim(1:3);
VO.mat        = mat;
VO            = spm_create_vol(VO);
T = [];
for i = 1:imgdim(3)
    M = inv(spm_matrix([0 0 -i])*inv(VO.mat)*V.mat);
    mask_slice = spm_slice_vol(V,M,VO.dim(1:2),1);
    ind0 = find(feval('gt', mask_slice,0));
    img = spm_slice_vol(V, M, imgdim(1:2), 1);
    if ~isempty(ind0), 
        T = [T; img(ind0) ]; 
    end
end
T=T';

