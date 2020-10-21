function [f,g] = ml_ordinal_logistic_loss(w,X,y,k)
% ml_ordinal_logistic_loss(w,X,y,k)
%
% Description:
%   - Calculate negative log likelihood (nll) and gradient w.r.t. weights
%     and threshold values for ordinal logistic regression.
%      
% Inputs:
%   - w(nFeatures,1): a weight vector
%   - X(nTrain, nFeatures): a matrix containing the training data
%   - y(nTrain, 1): a vector containing the target variable for each
%                   training example
%   - k: number of orders
%
% Returns:
%   - f: the nll function value 
%   - g: the gradient value of the nll 
%
% Authors:
%   - Mark Schmidt (2014); adapted for matLearn by Geoffrey Roeder (2016)

[n,p] = size(X);
nVars = length(w);
gamma = [-inf;0;cumsum(w(p+1:end));inf];
w = w(1:p);

% F(x) := sigmoid(x)
L = F(gamma(y+1) - X*w) - F(gamma(y) - X*w);
f = -sum(log(L));

if nargout > 1
   g = zeros(nVars,1);
   
   % Derivative wrt weights
   sigm1 = F(gamma(y+1) - X*w);
   sigm2 = F(gamma(y) - X*w);
   inner = (sigm1.*(1-sigm1) - sigm2.*(1-sigm2))./L;
   g(1:p) = X'*inner;

   % Derivative wrt cutoffs
   for i = 1:n
      if y(i) >= 3
         g(p+1:p+y(i)-2) = g(p+1:p+y(i)-2) + sigm2(i)*(1-sigm2(i))/L(i);
      end
      if y(i)+1 >= 3 && y(i)+1 <= k
          g(p+1:p+y(i)+1-2) = g(p+1:p+y(i)+1-2) - sigm1(i)*(1-sigm1(i))/L(i);
      end
   end
end
end

function [sigm] = F(x)
    sigm = 1./(1+exp(-x));
end