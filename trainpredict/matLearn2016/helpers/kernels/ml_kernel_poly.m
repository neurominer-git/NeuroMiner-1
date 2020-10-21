function [kern, name] = ml_kernel_poly(X1,X2,options)
% ml_kernel_poly(X1,X2,options)
%
% Description:
%   - Generates Gram matrix for feature vectors X1 and X2 in polynomial  
%     basis of specified order
%
% Options:
%   - bias: 1 to add a bias term to polynomial basis, 0 otherwise
%   - order: degree of polynomial basis


[bias, order] = myProcessOptions(options, 'bias', 0, 'order', 2);
kern = (X1*X2'+bias).^order;
name = 'Polynomial Kernel';
end