function [ fY, C ]= nk_SpatFilt(Y, SPATIAL)

if isfield(SPATIAL,'PX') && numel(SPATIAL.PX.opt)>1
    n_opt = numel(SPATIAL.PX.opt);
else
    n_opt=1;
end
fY = cell(n_opt,1); 

if ~isempty(SPATIAL.badcoords);
    tY = zeros(size(Y,1),size(SPATIAL.badcoords,2));
    tY(:,~SPATIAL.badcoords) = Y;
else
    tY = Y;
end

for i=1:n_opt
    switch SPATIAL.cubetype
        case 1
            fY = Y; C = []; return
        case 2
            fprintf('\nAbs. diff. weighting (6 neighbors) \n')
        case 3
            fprintf('\nSpat. var. weighting (27 neighbors) \n')
        case 4
            if isfield(SPATIAL,'cubefwhm'), fwhm  = SPATIAL.PX.opt(i); end
            fprintf('\nGaussian smoothing with %g FWHM\n', fwhm(1))
        case 5
            if isfield(SPATIAL,'cubevoxdim'), voxdim = SPATIAL.PX.opt(i); end
            fprintf('\nReslicing to %g mm isotropic voxels\n', voxdim(1))
    end 
    C = zeros(size(tY));
    for j=1:size(Y,1)
        switch SPATIAL.cubetype
            case 2
                
                V = zeros(SPATIAL.dims);
                V(SPATIAL.indvol) = tY(j,:); % transfer vector into 3D space

                [m,n,p] = size(V);
                rowC = 1:m; rowN = [1 1:m-1]; rowS = [2:m m];
                colC = 1:n; colE = [1 1:n-1]; colW = [2:n n];
                sliC = 1:p; sliD = [1 1:p-1]; sliU = [2:p p];

                % Compute difference between center pixel and each of the voxel 6 neighbors.
                north = V(rowN,colC,sliC)-V(rowC,colC,sliC);
                south = V(rowS,colC,sliC)-V(rowC,colC,sliC);
                east = V(rowC,colE,sliC)-V(rowC,colC,sliC);
                west = V(rowC,colW,sliC)-V(rowC,colC,sliC);
                up = V(rowC,colC,sliU)-V(rowC,colC,sliC);
                down = V(rowC,colC,sliD)-V(rowC,colC,sliC);
               
                % Compute sum of absolute differences
                T = abs(north) + abs(south) + abs(east) + abs(west) + abs(up) + abs(down);
                C(j,:) = T(SPATIAL.indvol);
               
            case 3

                V = zeros(SPATIAL.dims);
                V(SPATIAL.indvol) = tY(j,:); % transfer vector into 3D space
                T = spatconst27(V);
                C(j,:) = log(T(SPATIAL.indvol));
                C(j,:) = C(j,:)./max(C(j,:)); % scale to range [0,1]
               
            case 4
                 
                V = zeros(SPATIAL.dims);
                V(SPATIAL.indvol) = tY(j,:); % transfer vector into 3D space
                T = zeros(size(V));
                if fwhm
                    nk_smooth(SPATIAL.Vm, V, T, fwhm);
                else
                    T = V;
                end
                C(j,:) = T(SPATIAL.indvol); 
                
            case 5
                
                if voxdim
                    Vm = SPATIAL.Vm;
                    V = zeros(SPATIAL.dims);
                    V(SPATIAL.indvol) = tY(j,:); % transfer vector into 3D space
                    T = nk_ResizeImg(Vm, repmat(voxdim,1,3), nan(2,3));
                else
                    T = tY(j,:);
                end
                if j==1,
                    C = zeros(size(tY,1),size(T,2));
                end
                C(j,:) = T;
        end
        fprintf('.');
    end
    
    switch SPATIAL.cubetype
        case {2, 3}
            xY = tY./C;
        case {4, 5}
            xY = C;
    end
    
    if ~isempty(SPATIAL.badcoords) && SPATIAL.cubetype ~= 5; xY = xY(:,~SPATIAL.badcoords); end
    fY{i} = xY;
    
end