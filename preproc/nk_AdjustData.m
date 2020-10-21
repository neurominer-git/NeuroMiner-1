function nk_AdjustData

spm5ver = 0; if strcmp(spm('ver'),'SPM5') || strcmp(spm('ver'),'SPM8') , spm5ver = 1; end

if spm5ver
    brainmask = spm_select(1,'image','select brainmask');
else
    brainmask = spm_get(1,'IMAGE','select brainmask');
end

Vm  = spm_vol(brainmask);

if spm5ver
    P = spm_select(Inf,'image','Select images');
else	
    P = spm_get(Inf,'IMAGE','Select images');
end

if isempty(P), return; end;

% Load volume infos
V = spm_vol(P);
n = size(V,1);

G = nk_input('Nuisance variables',0,'r',[],[n Inf]);

% Setup for Proportional Scaling
GlobNorm = nk_input('Global normalisation ?',0,'y/n',[1,0],2);

if GlobNorm
    GM = nk_input('PropScaling to ',0,'e',100);
else
    GM = 1;
end

% Scale data to some global value
scaling = nk_input('Global scaling',0,'m','None|User specified globals|Compute as mean voxel value',[1 2 3],1);

switch scaling
    
    case 1
        if ~GM 
            for i = 1:n
                V(i).pinfo(1:2,:) = V(i).pinfo(1:2,:)*100/GM;
            end
        end
    case 2
        g = nk_input(['Globals '],0,'r',[],[n,1]);
        gSF    = GM./g;
     
        for i = 1:n
            V(i).pinfo(1:2,:) = V(i).pinfo(1:2,:)*gSF(i);
        end

    case 3	
        %-Compute as mean voxel value (within per image fullmean/8 mask)
        fprintf('Calculating globals\n')
        for i = 1:n
            g = spm_global(V(i));
            gSF = GM/g;
            V(i).pinfo(1:2,:) = V(i).pinfo(1:2,:)*gSF;
        end
        
end

suff = nk_input('Suffix for adjusted files',0,'s',['_adj']);

indvol=[];Y=[];
fprintf('\nConcatenating data into Y.')
for sl=1:Vm.dim(3) % Loop through z dimension
   
    % read mask
    M = spm_matrix([0 0 sl 0 0 0 1 1 1]);
    M1  = Vm.mat\V(1).mat\M;
    mask_slice = spm_slice_vol(Vm,M,Vm.dim(1:2),1);
    fprintf(1,'\nReading slice %g/%g: ',sl,Vm.dim(3))
    ind0 = find(mask_slice > 0.5);
    ind = ind0 + (sl - 1)*prod(Vm.dim(1:2));
    indvol = [indvol; ind];
    clear mask_slice

    % read data inside masked slice
    if ~isempty(ind0)
        y = zeros(n, length(ind0));
        for i = 1:n
            fprintf(1,'.')	
            d = spm_slice_vol(V(i),M1,Vm.dim(1:2),1);
            y(i,:) = d(ind0);
        end
        Y = [Y y];
    end
end
fprintf('\nFinished concatenating data into Y.')

% Remove Nuisance effects
Y = nk_PartialCorrelations(G,Y);

% Write adjusted data
for i=1:n
    
    [pth,nam,ext] = fileparts(deblank(P(i,:)));
    if strcmp(ext,'.nii,1')
        ext = regexprep(ext,'.nii,1','.nii');
    else
        ext = regexprep(ext,'.img,1','.img');
    end
    V(i).fname = fullfile(pth,[nam suff ext]);
    V(i).descrip = ['Adjusted ' V(i).descrip];
    Vimg = zeros(Vm.dim(1:3));
    Vimg(indvol) = Y(i,:);
    fprintf('\nWriting %s to disk.',V(i).fname)
    spm_write_vol(V(i),Vimg);
    
end
fprintf('\nDone.\n')