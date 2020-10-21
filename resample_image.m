voxsiz = [1 1 1]; % new voxel size {mm}
V = spm_select([1 Inf],'image');
V = spm_vol(V);
VV = [V(1) V(1)];
for i=1:numel(V)
   bb        = spm_get_bbox(V(i));
   VV(1:2)   = V(i);
   VV(1).mat = spm_matrix([bb(1,:) 0 0 0 voxsiz])*spm_matrix([-1 -1 -1]);
   VV(1).dim = ceil(VV(1).mat \ [bb(2,:) 1]' - 0.1)';
   VV(1).dim = VV(1).dim(1:3);
   spm_reslice(VV,struct('mean',false,'which',1,'interp',1)); % 1 for linear
end