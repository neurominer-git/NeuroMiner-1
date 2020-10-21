function [sY, VO] = nk_PerfResampleImages(Vm, Y, Voxdim, BB, interp)

nY = 1; if iscell(Y), nY = numel(Y); sY = cell(nY,1); end
if ~exist('Voxdim','var') || isempty(Voxdim), Voxdim =[nan nan nan]; end
if ~exist('BB','var') || isempty(BB), BB =[nan nan nan; nan nan nan]; end
if ~exist('interp','var') || isempty(interp), interp = 1; end

if any(isnan(Voxdim))
    vprm = spm_imatrix(Vm.mat);
    vvoxdim = vprm(7:9);
    Voxdim(isnan(Voxdim)) = vvoxdim(isnan(Voxdim));
end
Voxdim = Voxdim(:)';

mn = BB(1,:);
mx = BB(2,:);

% default BB to current volume's
if any(isnan(BB(:)))
    vbb = world_bb(Vm);
    vmn = vbb(1,:);
    vmx = vbb(2,:);
    mn(isnan(mn)) = vmn(isnan(mn));
    mx(isnan(mx)) = vmx(isnan(mx));
end

% voxel [1 1 1] of output should map to BB mn
% (the combination of matrices below first maps [1 1 1] to [0 0 0])
mat = spm_matrix([mn 0 0 0 Voxdim])*spm_matrix([-1 -1 -1]);
% voxel-coords of BB mx gives number of voxels required
% (round up if more than a tenth of a voxel over)
imgdim = ceil(mat \ [mx 1]' - 0.1)';

VO = [];
if nargout == 2 % Resampled output image ?
    VO            = Vm;
    [pth,nam,ext] = fileparts(Vm.fname);
    VO.fname      = fullfile(pth,['r' nam ext]);
    VO.dim(1:3)   = imgdim(1:3);
    VO.mat        = mat;
end

for j = 1:nY
    if nY>1, Yj = Y{j}; else Yj=Y; end
    img = zeros(imgdim(1:3));
    for i = 1:imgdim(3)
        M = inv(spm_matrix([0 0 -i])*inv(mat)*Vm.mat);
        img(:,:,i) = spm_slice_vol(Yj, M, imgdim(1:2), interp); % 
    end
    if nY>1, sY{j} = img; else sY = img; end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function bb = world_bb(V)
%  world-bb -- get bounding box in world (mm) coordinates

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

% bounding box (world) min and max
mn = min(tc,[],2)';
mx = max(tc,[],2)';
bb = [mn; mx];