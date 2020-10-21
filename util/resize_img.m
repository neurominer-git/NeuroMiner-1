function resize_img(imnames, Voxdim, BB)
%  resize_img -- resample images to have specified voxel dims and BBox
% resize_img(imnames, voxdim, bb)
%
% Output images will be prefixed with 'r', and will have voxel dimensions
% equal to voxdim. Use NaNs to determine voxdims from transformation matrix
% of input image(s).
% If bb == nan(2,3), bounding box will include entire original image
% Origin will move appropriately.

% Based on John Ashburner's reorient.m
% http://www.sph.umich.edu/~nichols/JohnsGems.html#Gem7
% http://www.sph.umich.edu/~nichols/JohnsGems5.html#Gem2

% Check spm version:
if exist('spm_select','file') % should be true for spm5
    spm5 = 1;
    select = @(msg) spm_select(inf, 'image', msg);
elseif exist('spm_get','file') % should be true for spm2
    spm5 = 0;
    select = @(msg) spm_get(inf, 'img', msg);
else
    error('Failed to locate spm_get or spm_select; please add SPM to Matlab path')
end

% set up defaults (including analyze.flip)
spm_defaults;

% prompt for missing arguments
if ( ~exist('imnames','var') || isempty(char(imnames)) )
    imnames = select('Choose images to resize');
end
% check if inter fig already open, don't close later if so...
Fint = spm_figure('FindWin', 'Interactive'); Fnew = [];
if ( ~exist('Voxdim', 'var') || isempty(Voxdim) )
    Fnew = spm_figure('GetWin', 'Interactive');
    Voxdim = spm_input('Vox Dims (NaN for "as input")? ', '+1', 'e', '[nan nan nan]', 3);
end
if ( ~exist('BB', 'var') || isempty(BB) )
    Fnew = spm_figure('GetWin', 'Interactive');
    BB = spm_input('Bound Box (NaN => original)? ', '+1', 'e', '[nan nan nan; nan nan nan]', [2 3]);
end

% reslice images one-by-one
vols = spm_vol(imnames);
for V=vols'
    % (copy to allow defaulting of NaNs differently for each volume)
    voxdim = Voxdim;
    bb = BB;
    % default voxdim to current volume's voxdim, (from mat parameters)
    if any(isnan(voxdim))
        vprm = spm_imatrix(V.mat);
        vvoxdim = vprm(7:9);
        voxdim(isnan(voxdim)) = vvoxdim(isnan(voxdim));
    end
    voxdim = abs(voxdim(:))'; % (handle analyze_flip separately)
    
    bbmn = bb(1,:);
    bbmx = bb(2,:);
    % default BB to current volume's
    if any(isnan(bb(:)))
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
    mat = spm_matrix([bbmn 0 0 0 voxdim])*spm_matrix([-1 -1 -1]);
    % voxel-coords of BB mx gives number of voxels required (now round up)
    imgdim = ceil(mat \ [bbmx 1]')';

    % reflect in x if required
    if spm_flip_analyze_images; mat = diag([-1 1 1 1])*mat; end;
    
    % output image
	VO            = V;
	[pth,nam,ext] = fileparts(V.fname);
	VO.fname      = fullfile(pth,['r' nam ext]);
	VO.dim(1:3)   = imgdim(1:3);
	VO.mat        = mat;
    VO = spm_create_vol(VO);
    spm_progress_bar('Init',imgdim(3),'reslicing...','planes completed');
    for i = 1:imgdim(3)
		M = inv(spm_matrix([0 0 -i])*inv(VO.mat)*V.mat);
		img = spm_slice_vol(V, M, imgdim(1:2), 1);
		spm_write_plane(VO, img, i);
		spm_progress_bar('Set', i)
    end
	spm_progress_bar('Clear');
end
if (isempty(Fint) && ~isempty(Fnew))
    % interactive figure was opened by this script, so close it again.
    close(Fnew);
end
disp('Done.')


