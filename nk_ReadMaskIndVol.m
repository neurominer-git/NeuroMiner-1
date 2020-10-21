function [dims, indvol, Y, vox] = nk_ReadMaskIndVol(Vm, Vp, label, labelop)
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% function [dims,indvol] = nk_ReadMaskIndVol(Vm)
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

indvol=[];
if ~exist('Vm','var'),
    Pm = spm_select(1,'image','Space-defining image');
    Vm = spm_vol(Pm);
elseif exist('Vm','var') &&  ~isstruct(Vm) && exist(Vm,'file')
    Vm = spm_vol(Vm);
end

if ~exist('Vp','var'),
    Pp = spm_select(1,'image','Image to read');
    Vp = spm_vol(Pp);
elseif exist('Vp','var') && ~isstruct(Vp) && exist(Vp,'file')
     Vp = spm_vol(Vp);
end

if ~exist('label','var') || isempty(label), label = 0; end
if ~exist('labelop','var') || isempty(labelop), labelop = 'gt'; end

vox = sqrt(sum(Vm.mat(1:3,1:3).^2));                
Y   =[];
for sl=1:Vm.dim(3)
    % read mask
   
    M = spm_matrix([0 0 sl 0 0 0 1 1 1]);
   
    mask_slice = spm_slice_vol(Vm,M,Vm.dim(1:2),1);
    ind0 = find(feval(labelop, mask_slice,label));
    ind = ind0 + (sl - 1)*prod(Vm.dim(1:2));
    if exist('Vp','var') && ~isempty(Vp)
        M1 = Vm.mat\Vp.mat\M;
        d = spm_slice_vol(Vp,M1,Vm.dim(1:2),1);
        if ~isempty(d(ind0))
            Y = [Y; d(ind0)];
        end
    end
    indvol = [indvol; ind];
    clear mask_slice
end

dims = Vm.dim(1:3);