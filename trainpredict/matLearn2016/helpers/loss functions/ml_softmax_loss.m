function [nll,g,H] = ml_softmax_loss(w,X,y,k)
% [nll,g,H] = ml_softmax_loss(w,X,y,k)
% 
% Description: Calculates softmax loss where weights for last class are 
%              fixed at 0 to avoid overparameterization
%
% Parameters:
%   - w(feature*class,1): weights (assumed to be 0 for last class)
%   - X(nInstances,nFeatures): observed features 
%   - y(nInstances,1): target labels
%   - k: number of classes


[n,d] = size(X);
w = reshape(w,[d k-1]);
w(:,k) = zeros(d,1);

Z = sum(exp(X*w),2);
nll = -sum((sum(X.*w(:,y).',2) - log(Z)));

if nargout > 1
    g = zeros(d,k-1);

    for c = 1:k-1
        g(:,c) = -sum(X.*repmat((y==c) - exp(X*w(:,c))./Z,[1 d]));
    end
    g = reshape(g,[d*(k-1) 1]);
end

if nargout > 2
    H = zeros(d*(k-1));
    SM = exp(X*w(:,1:k-1))./repmat(Z,[1 k-1]);
    for c1 = 1:k-1
        for c2 = 1:k-1
            D = SM(:,c1).*((c1==c2)-SM(:,c2));
            H((d*(c1-1)+1):d*c1,(d*(c2-1)+1):d*c2) = X'*diag(sparse(D))*X;
        end
    end
end
