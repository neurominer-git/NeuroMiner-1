function [model] = ml_regression_regressOnOne(X,y,options)
% ml_regression_regressOnOne(X,y,options)
%
% Description:
%    - Minimizes the squared error for one feature
%
% Options:
%    - selectedFeature: The feature to use for the univariate regression
%
% Authors:
%    - Mark Schmidt (2014)

[nTrain,nFeatures] = size(X);

[selectedFeature] = myProcessOptions(options,'selectedFeature',1);

% Compute the regression weight w for feature j, ignoring the others
j = 1;
x2 = 0;
xy = 0;
for i = 1:nTrain
   x2 = x2 + X(i,j)^2;
   xy = xy + X(i,j)*y(i);
end
w = xy/x2;

model.name = 'Regress on One';
model.selectedFeature = selectedFeature;
model.w = w;
model.predict = @predict;

end

function [yhat] = predict(model,Xhat)
    [nTest,nFeatures] = size(Xhat);
    w = model.w;
    j = model.selectedFeature;
    
    yhat = zeros(nTest,1);
    for i = 1:nTest
        yhat(i) = Xhat(i,j)*w;
    end
end
