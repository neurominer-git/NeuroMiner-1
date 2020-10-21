function [kern, name] = ml_kernel_rbf(X1,X2,options)

[sigma] = myProcessOptions(options, 'sigma', 1);
[n1, p] = size(X1);
n2 = size(X2,1);
Z = 1/sqrt(2*sigma^2);

% Vectorized way from Alex Smola's blog
D = X1.^2*ones(p,n2) + ones(n1,p)*(X2').^2 - 2*X1*X2';
kern = Z*exp(-D/(2*sigma^2));

name = 'RBF Kernel';
end