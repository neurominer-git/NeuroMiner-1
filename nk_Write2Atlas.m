function Vm = nk_Write2Atlas(Y, fname, Vm)

if ~exist('Vm','var') || ~exist(Vm, 'file')
    Vm = spm_select(1,'image','Select atlas image');
end
Vmvol = spm_vol(Vm);
A = spm_read_vols(Vmvol);

Al = unique(A(:));
if Al(1) == 0, Al(1)=[]; end
if numel(Al) ~= numel(Y)
    error('Vector Y does not have the same dimensionality as the labels in the atlas')
end
T = A;
for i=1:numel(Al)
   
    ind = A == i;
    T(ind) = Y(i);
   
end

fname = fullfile(pwd,[fname '.nii']);
Vmvol.fname = fname;
Vmvol.dt = [16,0];
spm_write_vol(Vmvol,T)