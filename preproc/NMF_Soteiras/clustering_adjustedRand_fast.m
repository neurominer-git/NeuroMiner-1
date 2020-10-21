function [r]=clustering_adjustedRand_fast(u,v)

% clustering quality measures assumptions : 
% 0 corresponds to background
% we do not care about the background
% cluster labels are assumed to be enumerated from 1 to max number of
% clusters
%
% this function should not be used when comparing binary segmentations

m=max(max(u),max(v));

if(m == 1)
    error('ClusteringAdjustedRandFast:argChk','This method should not be used for comparing binary segmentations');
end

va=zeros(1,m);
vb=zeros(1,m);
mat=zeros(m);

for i=1:m
    va(i) = sum(u==i) ;
    vb(i) = sum(v==i) ;
    hu = (u==i) ;
    for j=1:m
        hv = (v==j) ;
        mat(i,j) = sum(hu.*hv);
    end
end

ra=sum(sum(mat.*(mat-ones(m))))/2.0;

rb=va*(va-ones(1,m))'/2.0;

rc=vb*(vb-ones(1,m))'/2.0;

rn=length(u)*(length(u)-1)/2.0;

r=(ra-rb*rc/rn)/( 0.5*rb+0.5*rc-rb*rc/rn );


