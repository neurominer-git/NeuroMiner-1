function [model] = ml_generative_KDE(X,~,options)
% ml_generative_KDE(X,y,options)
% 
% Description:
%    - Computes the maximum likelihood kernel density estimate fit for the data. Targets y
%       are expected to be empty and are ignored. The prediction function
%       returns an [nTest, 1] vector with the likelihood of each data point in
%       the test set.
% 
% Options:
%    - kernelFunc: isotropic kernels for KDE (default: rbf) 
%    - kernelOptions: options for kernelFunc (default: sigma = 1)
%    - weights: vector to weight training data (default: vector of 1's)
% 
% Authors:
% 	 - Sharan Vaswani (2014): sharanv@cs.ubc.ca

[nTrain, nFeatures] = size(X);

kernel_options = [];
kernel_options.sigma = 1;
[kernelFunc,kernelOptions,z] = myProcessOptions(options,...
    'kernelFunc',@ml_kernel_rbf, ...
    'kernelOptions',kernel_options, ...
    'weights', ones(nTrain,1));

model.X = X;
model.kernelFunc = kernelFunc;
model.kernelOptions = kernelOptions;
model.weights = z;
model.name = 'Generative Kernel Density Estimation Model';
model.predict = @predict;
end

function [lik] = predict(model,Xhat)
% Prediction function
[nTrain,nFeatures] = size(model.X);
[nTest,nFeatures] = size(Xhat);

lik =  (model.kernelFunc(Xhat,model.X,model.kernelOptions) * model.weights)/sum(model.weights);
end
