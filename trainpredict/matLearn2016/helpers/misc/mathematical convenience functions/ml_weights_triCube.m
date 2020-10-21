function [weights] = ml_weights_triCube(distances)
% ml_weights_triCube(distances)
%
% Description:
%    - Calculate respective weights of each neighbours given their distances
%       to test point, using tricube weighting function
% Inputs:
%    - distances: a vector of size k (# of neighbors)
%       distances and return a vector of size k of weights. The distances given
%       will be between 0 and 1, where 0 represents the point itself and 1
%       represents the farthest point
%
% Authors:
%   - Giorgio Gori (2014)

weights = (1 - distances.^3).^3;
end