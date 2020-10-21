function [model] = ml_multilabel_mutlinomial(X,y,options)
% ml_multilabel_multinomial(X,y,options)
%
% Description:
%	 - Multilabel classification using multiclass logistic regression
%
%
% Options:
%   - nLabels: number of labels in dataset (default: None)
%   - verbose: 1 for detailed optimizer output, 0 no output (default: 1)
%   - classifier: function handle for classifier
%
% Authors:
% 	 - Mark Schmidt (2014); adapted for matLearn by Geoffrey Roeder (2016)

[n,p] = size(X);

if nargin < 3
    options = [];
end

[nLabels,verbose,classifier] = ...
    myProcessOptions(options,'nLabels',[],'verbose',1,'classifier',...
                     @ml_multiclass_logistic);

if isempty(nLabels)
    assert(1==0,'nLabels must be specified for multi-label classifiers\n');
end

model.nLabels = nLabels;
model.nClasses = 2^nLabels;
model.subModel = classifier(X,y,options);
model.predict = @predict;
model.lossFunc = @loss;
model.name = 'Multilabel Multinomial';

end

function [y] = predict(model,X)
y = model.subModel.predict(model.subModel,X);
end

function [nll] = loss(model,X,y)
    nll = model.subModel.lossFunc(model.subModel,X,y);
end
