function rotatemat = bootprocrust(origlv,bootlv)
%syntax rotatemat=rri_bootprocrust(origlv,bootlv)

%define coordinate space between orignal and bootstrap LVs
temp=origlv'*bootlv;

%orthogonalze space
[V,W,U]=svd(temp);

%determine procrustean transform
rotatemat=U*V';
