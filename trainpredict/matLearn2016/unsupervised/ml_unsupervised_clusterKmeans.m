function [model] = ml_unsupervised_clusterKmeans(X, options)
% ml_unsupervised_clusterKmeans(X,options)
%
% Description:
%	 - Cluster the observations in X using Kmeans algorithm
%
% Options:
%   - K: target number of clusters. (Required parameter) 
%   - means: means initialization matrix
%   - kpp: 0 for random initialization, 1 for Kmeans++ initialization. Has
%          no effect if an initial means vector is provided.
%   - plot: 1 for visualization plot, 0 for no adaptive plot. Default: 1.
%
% Authors:
% 	 - Mark Schmidt (2014)
%	 - Geoffrey Roeder (2016)


if nargin < 2
    y = NaN;
    options = [];
end

[nInstances,nVars] = size(X);
[k, means, kpp, verbose, plot] = myProcessOptions(options, 'k', -1, ...
                                                 'means', [], 'kpp', 1, ...
                                                 'verbose', 0, ...
                                                 'plot', 1);

if k == -1
    error('Must specify a target number number of clusters, k');
end

if isempty(means)
    if(kpp)
        % K-means++ improved initialization
        means = X(ceil(rand*nInstances),:);
        
        for k = 2:k
            % Compute squared distance of every point to the closest mean
            distances = X.^2*ones(nVars,k-1)...
                        + ones(nInstances,nVars)*(means').^2-2*X*means';
            minDistances = min(distances,[],2);
            
            % Sample the next proportional to these distances squared
            i = sampleDiscrete(minDistances/sum(minDistances));
            means(k,:) = X(i,:);
        end
        initStr = 'K-means++';
    else
        % Choose random points to initialize means
        nSamples = 1;
        for k = 1:k
            samples = ceil(rand(nSamples,1)*nInstances);
            means(k,:) = mean(X(samples,:),1);
        end
        initStr = 'random';
    end
else
    initStr = 'user-specified';
end

X2 = X.^2*ones(nVars,k);
while 1
    means_old = means;
    
    % Compute Euclidean distance between each data point and each mean
    D = X2 + ones(nInstances,nVars)*(means').^2 - 2*X*means'; 
    
    % Assign each data point to closest mean
    [~,assignment] = min(D,[],2);
    
    % Compute mean of each cluster
    means = zeros(k,nVars);
    for k = 1:k
        means(k,:) = mean(X(assignment==k,:),1);
    end
    
    if nVars == 2 && plot
        % Make plot
        clf;hold on;
        colors = getColors;
        for k = 1:k
            h = plot(X(assignment==k,1),X(assignment==k,2),'.');
            set(h,'Color',colors{k});
        end
        pause(.25);
    end
    if verbose
     fprintf('Running K-means, difference = %f\n', ...
             max(abs(means-means_old)));
    end
    
    if max(abs(means-means_old)) < 1e-5
        break;
    end
end

model.means = means;
model.assignment = assignment;
model.predict = @predict;
model.name = strcat(['K-means Clustering ( ',initStr,' initialization)']);
end

function [clusters] = predict(model,Xhat)
% Returns closest cluster for each new point
% NOTE: this function does not recompute cluster means
k = length(unique(model.gamma));
means = model.means;
[nTest,nVars] = size(Xhat);
Xhat2 = Xhat.^2*ones(nVars,k);
D = Xhat2 + ones(nTest,nVars)*(means').^2 - 2*Xhat*means';
% Assign each data point to closest mean
[~,clusters] = min(D,[],2);
end

