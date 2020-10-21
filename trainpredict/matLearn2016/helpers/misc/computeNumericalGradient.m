function grad = computeNumericalGradient(func, w)
% numgrad = computeNumericalGradient(J, theta)
% w: a vector of parameters for a function
% func: Calling y = func(theta) must return the function value at theta. 
  
% preallocated numerical gradient
grad = zeros(size(w));

eps = 1e-4;
nVar = size(w, 1);
I = eye(nVar);
for i = 1:nVar 
    grad(i) = (func(w + eps * I(:, i)) - func(w - eps * I(:, i))) / ...
                 (2 * eps);
end
end
