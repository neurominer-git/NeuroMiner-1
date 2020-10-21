function indvol = nk_ExtrIndVol(Vm, thr, throp)

indvol=[];
for sl=1:Vm.dim(3)
    M = spm_matrix([0 0 sl 0 0 0 1 1 1]);
    mask_slice = spm_slice_vol(Vm,M,Vm.dim(1:2),1);
    ind0 = find(feval(throp, mask_slice, thr));
    ind = ind0 + (sl - 1)*prod(Vm.dim(1:2));
    indvol = [indvol; ind];
end
