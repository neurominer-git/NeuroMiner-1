function Consistency = nk_SpatConst(cubetype, Y, dims, indvol, fwhm)

[ix,kbin]=size(Y);
Consistency = zeros(ix,kbin);

switch cubetype

    case 1 %
        
        fprintf('spat. const. 4x4 ... ')
        for i=1:kbin
            V = zeros(dims);
            V(indvol) = Y(:,i); % transfer vector into 3D space

            [m,n,p] = size(V);
            rowC = 1:m; rowN = [1 1:m-1]; rowS = [2:m m];
            colC = 1:n; colE = [1 1:n-1]; colW = [2:n n];
            sliC = 1:p; sliD = [1 1:p-1]; sliU = [2:p p];

            % Compute difference between center pixel and each of the 4
            % nearest neighbors.
            north = V(rowN,colC,sliC)-V(rowC,colC,sliC);
            south = V(rowS,colC,sliC)-V(rowC,colC,sliC);
            east = V(rowC,colE,sliC)-V(rowC,colC,sliC);
            west = V(rowC,colW,sliC)-V(rowC,colC,sliC);
            up = V(rowC,colC,sliU)-V(rowC,colC,sliC);
            down = V(rowC,colC,sliD)-V(rowC,colC,sliC);

            C = abs(north) + abs(south) + abs(east) + abs(west) + abs(up) + abs(down);
            Consistency(:,i) = C(indvol);
            Consistency(:,i) = Consistency(:,i)./max(Consistency(:,i)); % scale to range [0,1]
        end
        
    case 2

        fprintf('spat. const. 9x3 ... ')
        for i=1:kbin
            V = zeros(dims);
            V(indvol) = Y(:,i); % transfer vector into 3D space
            C = spatconst27(V);
            Consistency(:,i) = C(indvol);
            %disp(num2str(max(Consistency(:,i))));
            Consistency(:,i) = Consistency(:,i)./max(Consistency(:,i)); % scale to range [0,1]
        end
        
    case 3
        fprintf('spat. const. Gaussian smoothing %g mm ...',fwhm)
        for i=1:kbin
            V = zeros(dims);
            V(indvol) = Y(:,i); % transfer vector into 3D space
            C = zeros(size(V));
            spm_smooth(V,C,fwhm); C =1./C;
            Consistency(:,i) = C(indvol); 
            Consistency(:,i) = Consistency(:,i)./max(Consistency(:,i));
        end
       
end

