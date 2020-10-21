function [sortedDist,index] = getBestScores(dist,k)
% get the nearest K values from a distance vector

% sort if needed
if k>1
    [sortedDist,index] = sort(dist,'descend'); % sorted from the greatest to the lowest
    sortedDist = sortedDist(1:k); % Get the nearest k elements
    index = index(1:k);   % Get the corresponding index in dist
else
    [sortedDist,index] = max(dist);
end