function [Z] = ml_visualize_ISOMAP(X, options)
% ml_visualize_ISOMAP(X,options)
%
% Description:
%	 - Learn a manifold in dataset by applying Multi-Dimensional Scaling
%      to shortest path between all points in a weighted graph built from 
%      the original data points.
%
% Options:
%   - minVariance: the minimum amount of variance to be explained by the
%                   selected components (default: 1)
%   - maxComponents: the largest number of eigenvectors to be retained
%                    in the final basis, subject to minVariance constraint
%                   (default: nVars of input X)
%   - disconnected: indicate whether the graph is connected.
%                   Distances set to infinite by shortest path algorithm 
%                   will be changed to the maximium distance in graph 
%                   (the geodesic distance)
%
% Authors:
%    - Geoffrey Roeder (2016)
% 	 - Mark Schmidt (2014)

[nInstances,nVars] = size(X);

[K,names,disconnected,optOptions,dim] = myProcessOptions(options,...
                                                         'K',2,...
                                                         'names',{},...
                                                         'disconnected',1,...
                                                         'optOptions',{},...
                                                         'dim',2);

D = ml_geodesic_distance(X, options);

% Initialize low-dimensional representation with PCA
[U, S, V] = svd(X);
W = V(:,1:dim)';
Z = X*W';

% check computability of gradient
[~,g] = stress(Z(:),D,names);
if sum(isinf(g)) > 0
    fprintf(['FAILED: couldn''t evaluate gradient. Try setting ', ...
            '''disconnected'' to 1.\n']);
    return
end

Z(:) = minFunc(@stress,Z(:),optOptions,D,names);
end

function [f,g] = stress(Z,D,names)

n = length(D);
k = numel(Z)/n;

Z = reshape(Z,[n k]);

f = 0;
g = zeros(n,k);
for i = 1:n
    for j = i+1:n
        % Objective Function
        Dz = norm(Z(i,:)-Z(j,:));
        s = D(i,j) - Dz;
        f = f + (1/2)*s^2;
        
        % Gradient
        df = s;
        dgi = (Z(i,:)-Z(j,:))/Dz;
        dgj = (Z(j,:)-Z(i,:))/Dz;
        g(i,:) = g(i,:) - df*dgi;
        g(j,:) = g(j,:) - df*dgj;
    end
end
g = g(:);

% Make plot if using 2D representation
if k == 2
    figure(3);
    clf;
    plot(Z(:,1),Z(:,2),'.');
    hold on;
    if ~isempty(names)
        for i = 1:n
            text(Z(i,1),Z(i,2),names(i,:));
        end
    end
end
end
