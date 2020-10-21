function [RotMat, rotatemat] = nk_Procrustes(TempMat,SourceMat)

% Define coordinate space between orignal and rotated components
temp=TempMat'*SourceMat;

% Orthogonalze space
[V,~,U]=svd(temp);

% Determine procrustean transform
rotatemat=U*V';

% Finally transform SourceMat to minimize difference to TempMat
RotMat = SourceMat * rotatemat;
