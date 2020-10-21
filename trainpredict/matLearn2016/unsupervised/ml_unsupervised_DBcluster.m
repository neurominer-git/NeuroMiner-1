function [model] = ml_unsupervised_DBcluster(X, options)
% ml_unsupervised_DBcluster(X,options)
%
% Description:
%	 - Cluster the observations in X using DBSCAN where a cluster is
%      defined as a connected dense region of points
%
% Options:
%   - eps: minimum distance between points to be considered 'close' 
%          (default: 1) 
%   - minPts: number of 'close' points needed to define a cluster 
%             (default: 3)
%
% Authors:
% 	 - Mark Schmidt (2014)

if nargin < 2
    options = [];
end

[nInstances,nVars] = size(X);

[eps,minPts] = myProcessOptions(options,'eps',1,'minPts',3);

for i = 1:nInstances
    dist(i,i) = inf;
    for j = i+1:nInstances
        dist(i,j) = sum((X(i,:)-X(j,:)).^2);
    end
end
dist = dist+dist';

visited = zeros(nInstances,1);
gamma = zeros(nInstances,1);
K = 0;
for i = 1:nInstances
    if ~visited(i)
        visited(i) = 1;
        neighbors = find(dist(:,i) <= eps);
        if length(neighbors) >= minPts
            K = K + 1;
            [visited,gamma] = expand(X,i,neighbors,K,eps,minPts,dist,visited,gamma);
        end
        
        if nVars == 2
            % Make plot
            clf;hold on;
            colors = getColors;
            h = plot(X(gamma==0,1),X(gamma==0,2),'.');
            set(h,'Color',[0 0 0]);
            for k = 1:K
                h = plot(X(gamma==k,1),X(gamma==k,2),'.');
                set(h,'Color',colors{k});
            end
            pause(.25);
        end
    end
    
    
end
model.Xtrain = X;
model.gamma = gamma;
end

function [visited,gamma] = expand(X,i,neighbors,K,eps,minPts,dist,visited,gamma)
gamma(i) = K;
ind = 0;
while 1
    ind = ind+1;
    if ind > length(neighbors)
        break;
    end
    n = neighbors(ind);
    gamma(n) = K;
    
    if ~visited(n)
        visited(n) = 1;
        neighbors2 = find(dist(:,n) <= eps);
        if length(neighbors2) >= minPts
            neighbors = [neighbors;setdiff(neighbors2,neighbors)];
        end
    end
end
end