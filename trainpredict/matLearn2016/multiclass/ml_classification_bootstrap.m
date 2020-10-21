function [model] = ml_classification_bootstrap(X,y,options)
% w = ml_classification_bootstrap(X,y,k)
%
% k: number of classes
%
% Computes Bootstrapped Classifier

if nargin < 3
    options = [];
end

[nBootstraps,bootstrapClassifier, nClasses] = myProcessOptions(options, ...
    'nBootstraps',20,'bootstrapClassifier',@ml_generative_Gaussian, ...
    'nClasses',2);

[nInstances,nFeatures] = size(X);

model.nBootstraps = nBootstraps;
model.bootstrapClassifier = bootstrapClassifier;

for k = 1:nBootstraps
    if mod(k,100) == 0
        fprintf('Bootstrap %d of %d\n',k,nBootstraps);
    end
    
    sample = ceil(nInstances*rand(nInstances,1));
    
    model.subModel{k} = bootstrapClassifier(X(sample,:),y(sample),options);
end
model.nClasses = nClasses;
model.name = strcat(['Bootstrapped ', model.subModel{1}.name]);
model.predict = @predict;

end

function [y] = predict(model,X)

for k = 1:length(model.subModel)
   y(:,k) = model.subModel{k}.predict(model.subModel{k},X);
end
y = mode(y,2);
end