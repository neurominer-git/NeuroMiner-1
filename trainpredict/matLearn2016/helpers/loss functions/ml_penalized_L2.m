function [nll,g,H] = ml_penalized_L2(w,gradFunc,lambda,varargin)
% ml_penalized_L2(w,gradFunc,lambda,varargin)
%
% Description:
%   - Adds L2 regularization to loss/gradient function. Pass in extra
%     parameters for loss/gradient function as additional parameters after
%     lambda (varargin)
% Authors:
%   - Mark Schmidt (2014)

if nargout <= 1
    [nll] = gradFunc(w,varargin{:});
elseif nargout == 2
    [nll,g] = gradFunc(w,varargin{:});
else
    [nll,g,H] = gradFunc(w,varargin{:});
end

nll = nll+sum(lambda.*(w.^2));

if nargout > 1
    g = g + 2*lambda.*w;
end

if nargout > 2
    if isscalar(lambda)
        H = H + 2*lambda*eye(length(w));
    else
        H = H + diag(2*lambda);
    end
end