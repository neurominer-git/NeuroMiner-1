function [kern, name] = ml_kernel_gram(X1,X2,options)

kern = (X1*X2');
name = 'Linear Kernel';
end