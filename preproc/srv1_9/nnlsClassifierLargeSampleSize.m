function [testClassPredicted,sparsity,Y]=nnlsClassifierLargeSampleSize(trainSet,trainClass,testSet,testClass,option)
% NNLS Classifier: testSet=trainSet*Y, s.t. Y>=0.
% Usage:
% [testClassPredicted,sparsity]=nnlsClassifier(trainSet,trainClass,[],testClass)
% [testClassPredicted,sparsity]=nnlsClassifier(trainSet,trainClass,testSet,testClass)
% [testClassPredicted,sparsity]=nnlsClassifier(trainSet,trainClass,testSet,testClass,option)
% trainSet, matrix, the training set with samples in columns and features in rows.
% trainClass: column vector of numbers or string, the class labels of the traning set.
% testSet: matrix, the test set.
% testClass: column vector of numbers or string, the class labels of the
% test/unknown set. It is actually unused in this function, thus, set it [].
% option: struct, the options to configue this function:
% option.method, string, the optimization algorithm used to solve the NNLS problem. It could be
%     'nnls': used the NNLS algorithm (default);
%     'seminmfupdaterule': use the update rules based algorithm;
%     'sparsennls': used NNLS algorithm with sparse constraint.
% option.predicter: the method to find the class label of a test sample according to Y. It could be
%     'max': the same class label with the training sample with the maximum coefficient (default);
%     'kvote': select k training samples with the k largest coefficients, and decide the class labels by majority voting.
% option.k: scalar, only for option.predicter='kvote'. The default is 1.
% option.kernel, string, specifies the kernel. can be 'linear'(default),'polynomial','rbf','sigmoid','ds'
% option.param, scalar or column vector, the parameters for kernels, the default is [].
% testClassPredicted: column vector, the predicted class labels of the test/unknown samples.
% sparsity: scalar, the sparsity of the coefficient matrix Y.
% Contact Information:
% Yifeng Li
% University of Windsor
% li11112c@uwindsor.ca; yifeng.li.cn@gmail.com
% May 23, 2011

optionDefault.method='nnls';
optionDefault.predicter='subspace';
optionDefault.k=1;
optionDefault.kernel='linear';
optionDefault.param=[];
optionDefault.sampleSelMethod='knn'; % knn
optionDefault.kfold=10; % how many parts does the training data is divided
optionDefault.selectThreshold=0.7;
optionDefault.neighbor=5;
% optionDefault.search=false;
optionDefault.sparsityThreshold=1e-4;
if nargin<5
    option=optionDefault;
else
    option=mergeOption(option,optionDefault);
end
% if tensor data
if size(trainSet,3)>1
    trainSet=matrizicing(trainSet,3);
    testSet=matrizicing(testSet,3);
    trainSet=trainSet';
    testSet=testSet';
end

% % normalization to length 1
% trainSet=normc(trainSet);
% testSet=normc(testSet);

switch option.sampleSelMethod
    case 'nnls'
        [sampleSelected,classesSelected]=sampleSelNNLS(trainSet,trainClass,testSet,option);
    case 'knn'
        [sampleSelected,classesSelected]=sampleSelKNN(trainSet,trainClass,testSet,option);
end
fprintf('running nnls after sample selection...\n');
[testClassPredicted,sparsity,Y]=nnlsClassifier(sampleSelected,classesSelected,testSet,testClass,option);  
fprintf('nnls finished\n');
end

