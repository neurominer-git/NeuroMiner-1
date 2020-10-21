function [z] = ml_visualize_Sammon(X, options)

% ml_visualize_Sammon(X,options)
%
% Description:
%	 - Reduce dimensionality of X by minimizing weighted squared difference
%       in distances between pairs of points in original and reduced 
%       dimensional spaces. 
%
% Options:
%   - nComponents: the number of components for vectors in the lower
%                  dimensional space
%
% Authors:
% 	 - Mark Schmidt (2014)

if nargin < 2
    options = [];
end

[nInstances,nVars] = size(X);


[nComp] = myProcessOptions(options,'nComponents',nVars);

% Standardize Columns
[X,mu,sigma_data] = standardizeCols(X);

% Form distance matrix
D = zeros(nInstances);
for i = 1:nInstances
    for j = 1:nInstances
        D(i,j) = norm(X(i,:)-X(j,:));
    end
end

% Form initial low-dimensional representation
z = 1e-16*randn(nInstances,nComp);
%z = randn(nInstances,nComp);
options = [];
options.display = 0;
options.verboseI = 0;
z(:) = minFunc(@stress, z(:), options, D);

end

function [S,g] = stress(z,D)
visualize = 1;
nInst = length(D);
nComp = numel(z)/nInst;

z = reshape(z,nInst,nComp);

if nargout > 1
    g = zeros(size(z));
end

S = 0;
for i = 1:nInst
    for j = i+1:nInst
        nrmDist = norm(z(i,:)-z(j,:));
        s = D(i,j) - nrmDist;
        S = S + s^2/D(i,j);
        
        if nargout > 1
            g(i,:) = g(i,:) - 2*(s/nrmDist)*(z(i,:)-z(j,:))/D(i,j);
            g(j,:) = g(j,:) - 2*(s/nrmDist)*(z(j,:)-z(i,:))/D(i,j);
        end
    end
end

if nargout > 1
    g = g(:);
end

if visualize
    if nComp == 2
        figure(1);clf;
        plot(z(:,1),z(:,2),'.');
        pause(.01)
    else
        figure(2);clf;
        plot3(z(:,1),z(:,2),z(:,3),'.');
        pause(.01)
    end
end
end