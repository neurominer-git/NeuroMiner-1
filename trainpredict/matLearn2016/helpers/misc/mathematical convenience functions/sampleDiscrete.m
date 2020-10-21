function [y] = sampleDiscrete(p,r)
% Returns a sample from a discrete probability mass function indexed by p
% (assumes that p is already normalized)
if nargin < 2
    r = rand;
end
y = find(cumsum(p) > r,1);