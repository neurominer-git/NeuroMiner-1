function [model] = ml_binaryclass_brokenStump(X,y,options)
% ml_binaryclass_brokenStump(X,y,options)
%
% Description:
%    - Finds the best threshold across all features by minimizing the 
%      number of misclassifications
%
% Options:
%    - None
%
% Authors:
% 	 - Mark Schmidt (2014)

[nTrain,nFeatures] = size(X);

% Find best threshold
minErr = inf;
minVar = 0;
for j = 1:nFeatures
    thresholds = [sort(unique(X(:,j)));max(X(:,j))+eps];
    
    for t = thresholds'
        err = sum(X(y==-1,j) < t) + sum(X(y==1,j) >= t);
        
        if err < minErr
            minErr = err;
            selectedFeature = j;
            minThreshold = t;
        end
    end
end

% Model outputs
model.name = 'Broken Stump Binary Classification';
model.selectedFeature = selectedFeature;
model.threshold = minThreshold;
model.predict = @predict;
end

function yhat = predict(model,Xhat)
% Prediction function

index = Xhat(:,model.selectedFeature) < model.threshold;
yhat = sign(index - 0.5);

end
