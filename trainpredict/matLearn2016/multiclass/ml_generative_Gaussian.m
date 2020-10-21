function [model] = ml_generative_Gaussian(X,~,~)
% ml_multiclass_generativeGaussian(X,y,options)
% 
% Description:
%    - Computes the maximum likelihood Gaussian fit for the data. Targets y
%       are expected to be empty and are ignored. The prediction function
%       returns an [nTest, 1] vector with the likelihood of each data point in
%       the test set.
% 
% Options:
% 	 - None
% 
% Authors:
%  	 - Mark Schmidt
%    - Neil Traft (2014)

[nTrain, nFeatures] = size(X);

mu = mean(X,1)';
% Modify to be positive definite
sigma = cov(X)+1e-8*eye(nFeatures);

% Model ouputs
model.mu = mu;
model.sigma = sigma;
model.predict = @predict;
model.name = 'Generative Gaussian Model';
end

function [lik] = predict(model,X)
% Prediction function
mu = model.mu;
sigma = model.sigma;

% Form Gaussian distribution
lik = mvnpdf(X,mu',sigma);
end
