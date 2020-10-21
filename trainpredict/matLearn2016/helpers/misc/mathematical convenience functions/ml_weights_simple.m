function [weights] = ml_weights_simple(distances)
% ml_weights_simple(distances)
%
% Description:
%    - Calculate respective simple weights of each neighbours given
%       their distances to test point
%
% Inputs:
%    - distances: a vector of size k (# of neighbof
%       distances and return a vector of size k of weights. The distances given
%       will be between 0 and 1, where 0 represents the point itself and 1
%       represents the farthest point
%
% Authors:
%   - Giorgio Gori (2014)

weights = 1 - distances;
end