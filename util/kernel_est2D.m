function [N,Z]=kernel_est2D(x,y,z,X,Y,sigma)
% function [Z,N]=kernel_est2D(x,y,z,X,Y,sigma);
% Gives a kernel estimate of the underlying density and undlerlying function z 
% INPUT: 
%   x: x-value of observation (Px1 vector)
%   y: y-value of observation (Px1 vector)
%   z: z-value of observation (Px1 vector): if z is empty, only density is returned
%   X: Meshgrid representation of X-coordinates 
%   Y: Meshgrid representation of Y-coordinates 
%   sigma: [SDx SDy], standard deviation of gaussian kernel in units of x
%   and y
% OUTPUT:
%   Z Kernel estimate of unlerlying values Z
%   N density estimate of the process (x,y)
if (~isempty(z))
    [num_nodes,num_cols]=size(z);
else
    [num_nodes,num_cols]=size(x);
end;    

Kernel=nm_normpdf([0:0.05:3.5],0,1)';
for i=1:length(X(:))
    d=[(x-X(i))/sigma(1) (y-Y(i))/sigma(2)];
    dn=sqrt(sum(d.^2,2));
    indx=find(dn<=3.5);
    din=floor(dn(indx)./0.05)+1;
    w=Kernel(repmat(din,1,num_cols));
    N(i)=sum(w);
    if (~isempty(indx) & ~isempty(z))
        Z(i)=sum(z(indx,:).*w)./sum(w);
    else 
        Z(i)=NaN;
    end;
end;
if (~isempty(z))
    Z=reshape(Z,size(X));
end;
N=reshape(N,size(X));
