function C = nk_CalcSpatConst(Y,Vm)

indvol=[];

% Read brainmask into vector format
for sl=1:Vm.dim(3)
    % read mask
    M = spm_matrix([0 0 sl 0 0 0 1 1 1]);
    %M1  = Vm.mat\V.mat\M;
    mask_slice = spm_slice_vol(Vm,M,Vm.dim(1:2),1);
    ind0 = find(mask_slice > 0.5);
    ind = ind0 + (sl - 1)*prod(Vm.dim(1:2));
    indvol = [indvol; ind];
    clear mask_slice
end

[n,k] = size(Y);

C = zeros(Vm.dim(1:3));
m=27; fm = zeros(1,m);
sz = size(C);

% Build feature matrix of 27 neighboring voxels of j training samples
for z = 2:sz(3)-1
	for y=2:sz(2)-1
		for x=2:sz(1)-1
    		% Initiate feature matrix
			% Extract neighboring voxels to columns
			for j=1:n
				fm(j,:) = reshape(V{j}(x-1:x+1,y-1:y+1,z-1:z+1),1,27);
			end
			% Now compute spatial consistency
			sfm = sum(sum(fm,2));
            if ~sfm, C(x,y,z) = 1; continue;end
			r = sum(fm,2)./m;
			g = sfm/(m*n);
			c = sum(fm)./n;	
			MSR = (m*sum(r-g)^2)/(n-1);
			MSC = (n*sum(c-g))/(m-1);
			MSE = (sfm^2-(n-1)*MSR-(m-1)*MSC-m*n*g^2)/((m-1)*(n-1));
			C(x,y,z) = (MSR-MSE)/(MSR+(m-1)*MSE+m/n*(MSC-MSE));
			fprintf('\nC(x=%g,y=%g,z=%g)=%1.3f',x,y,z,C(x,y,z))
		end
	end
end