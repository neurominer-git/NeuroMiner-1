function [V, indvol, vox] = nk_WriteVol(Vdat, fname, dim, spacimg, badcoords, label, labelop, indvol)
% ================================================================================
% FORMAT nk_WriteVol(Vdat, fname, dim, spacimg, badcoords, label, labelop, indvol)
% ================================================================================

sz = size(Vdat);

if sz(dim) > 1 % more than one volume is available
   k=sz(dim);
else
   k=1;
end
if exist('spacimg','var') && ~isempty(spacimg)
   spacimg = regexprep(spacimg,',1','');
else
   spacimg = []; 
end

if isempty(spacimg) || ~exist(spacimg,'file')
    
    spm5ver = nk_CheckSPMver;

    if spm5ver
        Vm = spm_vol(spm_select(1,'image','Space-defining image'));
    else
        Vm = spm_vol(spm_get(1,'IMAGE','Space-defining image'));
    end
else
    Vm = spm_vol(spacimg);
end
if ~exist('label','var') || isempty(label), label = 0.5; end
if ~exist('labelop','var') || isempty(labelop), labelop = 'gt'; end

if ~exist('indvol','var') || isempty(indvol)
    
    indvol=[];

    % Read brainmask into vector format
    for sl=1:Vm.dim(3)
        % read mask
        M = spm_matrix([0 0 sl 0 0 0 1 1 1]);
        %M1  = Vm.mat\V.mat\M;
        mask_slice = spm_slice_vol(Vm,M,Vm.dim(1:2),1);
        ind0 = find(feval(labelop,mask_slice,label));
        ind = ind0 + (sl - 1)*prod(Vm.dim(1:2));
        indvol = [indvol; ind];
        clear mask_slice
    end

    if exist('badcoords','var') && ~isempty(badcoords), 
        indvol = indvol(~badcoords); 
    end
end

if k>1
    V = cell(k,1); 
end
vox = sqrt(sum(Vm.mat(1:3,1:3).^2));

z=k;

while z>0
   
    Vol             = Vm;
    Vol.dim         = Vm.dim(1:3);
    Vol.dt          = [64,0];
%     if dim==1
%         Vol.pinfo(1)    = (max(Vdat(k,:)) - min(Vdat(k,:)))/32768;
%     else
%         Vol.pinfo(1)    = (max(Vdat(:,k)) - min(Vdat(:,k)))/32768;
%     end
    Vimg            = zeros(Vm.dim(1:3));    
    sz=size(Vdat);
    if isequal(sz,Vol.dim(1:3))
        Vimg=Vdat;
    else
        if dim == 1 
            Vimg(indvol) = Vdat(z,:);
        elseif dim == 2
            Vimg(indvol) = Vdat(:,z);
        end
    end
    
    if ~isempty(fname)
        if iscell(fname)
            F = fname{z};
        else
            F = deblank(fname(z,:));
        end
        fprintf('\nWriting %s to disk...', F);
        Vol.fname   = [pwd filesep F '.nii'];
        Vol.descrip = 'Generated by nk_WriteVol';
        %fprintf('\nWriting %s to disk: ',F)
        spm_write_vol(Vol,Vimg);
    else
        if k>1
            V{z} = Vimg;
        else
            V = Vimg;
        end
    end
    fprintf('\tDone.');
    z=z-1;
    
end

