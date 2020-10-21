function [model] = ml_regression_mean(X,y,options)
% ml_regression_mean(X,y,options)
%
% Description:
%	 - Always predicts the mean of the y values
%
% Options:
%	 - None
%
% Authors: Mark Schmidt (2014)

[nTrain,nFeaturs] = size(X);

% Compute the mean;
meanValue = 0;
for i = 1:nTrain
    meanValue = meanValue + y(i)/nTrain;
end

model.name = 'Mean';
model.meanValue = meanValue;
model.predict = @predict;

end

function [yhat] = predict(model,Xhat)
    [nTest,nFeatures] = size(Xhat);

    yhat = zeros(nTest,1);
    for i = 1:nTest
        yhat(i) = model.meanValue;
    end
end