function [model] = ml_multilabel_independent(X,y,options)
% ml_multilabel_independent(X,y,options)
%
% Description:
%	 - Naive multilabel classification using independent binary classifiers
%          for each label
%
% Options:
%   - nLabels: number of labels in dataset (default: None)
%   - verbose: 1 for detailed optimizer output, 0 no output (default: 1)
%   - binaryClassifier: function handle for binary classifier function
%                       (default: ml_binaryclass_logistic)
%
% Authors:
% 	 - Mark Schmidt (2014); adapted for matLearn by Geoffrey Roeder (2016)

[n,p] = size(X);

if nargin < 3
    options = [];
end

[nLabels,verbose,binaryClassifier] = myProcessOptions(options,'nLabels',...
        [],'verbose',1,'binaryClassifier',@ml_binaryclass_logistic);

if isempty(nLabels)
    assert(1==0,'nLabels must be specified for multi-label classifiers\n');
end

y = linearInd2Binary(y,nLabels);

model.nLabels = nLabels;
for k = 1:nLabels
   model.subModel{k} = binaryClassifier(X,y(:,k),options);
end

model.nClasses = 2^nLabels; % number of classes combinations (for plotting)
model.predict = @predict;
model.lossFunc = @loss;
model.name = 'Independent Logistic Classifiers';

end

function [y] = predict(model,X)
for k = 1:model.nLabels
   y(:,k) = model.subModel{k}.predict(model.subModel{k},X); 
end
y = binary2LinearInd(y);
end

function [nll] = loss(model,X,y)
y = linearInd2Binary(y,model.nLabels);
nll = 0;
for k = 1:model.nLabels
    nll = nll + model.subModel{k}.lossFunc(model.subModel{k},X,y(:,k));
end
end
