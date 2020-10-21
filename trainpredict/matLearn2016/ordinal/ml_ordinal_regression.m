function [model] = ml_ordinal_regression(X,y,options)
% ml_ordinal_regression(X,y,options)
%
% Description:
%    - Regression using ordinal regression model
%
% Options:
%   - regressionFunc: function handle for regression function
%                 (default: @ml_regression_leastSquares)
%
% Authors:
% 	- Mark Schmidt (2014)

[regressionFunc, subOptions] = myProcessOptions(options,'regressionFunc', ...
                                    @ml_regression_leastSquares, ...
                                    'subOptions', {});
model.nClasses = options.nClasses;
model.subModel = regressionFunc(X,y,subOptions);
model.predict = @predict;
model.name = strcat(['Ordinal',' ',model.subModel.name]);
 
function [y] = predict(model,X)
k = model.nClasses;
y = model.subModel.predict(model.subModel,X);
y = round(min(k,max(y,1)));

