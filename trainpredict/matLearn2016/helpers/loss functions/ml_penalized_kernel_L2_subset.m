function [nll,g,H] = ml_penalized_kernel_L2_subset(w,K,subset,gradFunc,lambda,varargin)
% ml_ordinal_logistic_loss(w,X,y,k)
%
% Description:
%   - Calculate negative log likelihood (nll) and gradient value for
%     ordinal logistic regression (you can use this instead of always 
%     adding it to the loss function code). This version applies the
%     regularization to a subset of the variables
%      
% Inputs:
%   - w(nFeatures,1): a weight vector of size nFeatures
%   - K: 
%   - subset: a vector containing indices of the weights to regularize
%   - gradFunc: number of orders
%   - lambda: scalar value for L2 penalty
%   - varargin: the other arguments needed for gradFunc
%
% Returns:
%   - f: the penalized nll function value 
%   - g: the gradient of f
%   - H: the Hessian of f
%
% Authors:
%   - Mark Schmidt; adapted for matLearn by Geoffrey Roeder (2016)

if nargout <= 1
    [nll] = gradFunc(w,varargin{:});
elseif nargout == 2
    [nll,g] = gradFunc(w,varargin{:});
else
    [nll,g,H] = gradFunc(w,varargin{:});
end

nll = nll+sum(lambda*w(subset)'*K*w(subset));

if nargout > 1
    g(subset) = g(subset) + 2*lambda*K*w(subset);
end

if nargout > 2
    H(subset,subset) = H(subset,subset) + 2*lambda*K;
end