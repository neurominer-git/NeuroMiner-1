function [testClassPredicteds,classPerforms,conMats,tElapseds,OtherOutputs]=multiClassifiers(trainSet,trainClass,testSet,testClass,methods,options)
% trainSet, matrix, the training set with samples in columns and features in rows.
% Usuage:
% [testClassPredicteds,classPerforms,tElapseds,OtherOutputs]=multiClassifiers(trainSet,trainClass,testSet,testClass,methods)
% [testClassPredicteds,classPerforms,tElapseds,OtherOutputs]=multiClassifiers(trainSet,trainClass,testSet,testClass,methods,options)
% trainClass: column vector of numbers or string, the class labels of the traning set.
% testSet: matrix, the test set.
% testClass: column vector of numbers or string, the class labels of the
%       test/unknown set. It is actually unused in this function, thus, set it [].
% methods: column vector of cell strings, designate which one or more classifier to use. It could be
%     'ksrsc': kernel sparse coding based classification methods, including NNLS,l1NNLS, and l1LS;
%     'subdic': sub-dictionary learning methods;
%     'svmm': my implementation of SVM classifiers including C-SVM, nu-SVM, and l_1-norm SVM;
%     'hdlm': high dimensional linear machine (the hard margin l_2-regularized linear model);
%     'nmf', 'rnmf': the local learning methods / clustering-and-classification methods;
%     'svdd': SVDD based multi-class classification, the nearest neighbor methods;
%     'hkm': hierarchial model for multi-class classification;
%     'nc': kernel nearest centroid method;
%     'lrc': linear regression classifier;
%     'logistic': logistic regression;
%     'bayesian': Bayesian classifier;
%     'naiveBayes': Naive Bayes classifier;
%     'knn': K-NN classifier;
%     'lsvm': local SVM classifier;
%%%     'kernelfda': (for future extension);
%%%    'psdkernelfda': (for future extension).
%      For example, methods={nnls;nmf;svm};
% options: column vector of cells. options{i} is the option for methods{i}. Type "help classification" for more information.
% testClassPredicteds: matrix, testClassPredicteds(:,i) are the predicted class labels by methods{i}.
% classPerforms: matrix, classPerforms(i,:) is the class performance by methods{i}.
% tElapseds: column vector; tElapseds(i) is the computing time of methods{i}.
% OtherOutputs: column vector of cell. OtherOutputs{i} is the other otherOutput of methods{i}.
% See also "classfication".
% Contact Information:
% Yifeng Li
% University of Windsor
% li11112c@uwindsor.ca; yifeng.li.cn@gmail.com
% May 18, 2011

%  % if tensor data
% iftensor=false;
% if size(trainSet,3)>1
%     iftensor=true;
%     numR=size(trainSet,1);
%     numC=size(trainSet,2);
%     trainSet=matrizicing(trainSet,3);
%     testSet=matrizicing(testSet,3);
%     trainSet=trainSet';
%     testSet=testSet';
%  end
numClasses=numel(unique([trainClass;testClass]));
numberPerfMeasure=2+numClasses;
% if numClasses<=4
%     numberPerfMeasure=6;
% else
%     numberPerfMeasure=2+numClasses;
% end
testClassPredicteds=nan(numel(testClass),numel(methods));
classPerforms=nan(numel(methods),numberPerfMeasure);
conMats=nan(numClasses,numClasses,numel(methods));
OtherOutputs=cell(numel(methods),1);
tElapseds=nan(numel(methods),1);
for m=1:numel(methods)
    method=methods{m};
    if isempty(options{m})
        option=[];
    else
        option=options{m};
    end
%     if iftensor
%        option.numR=numR;
%        option.numC=numC;
%     end
    fprintf('running method: %s ......\n',method);
    [model,OtherOutputTr]=classificationTrain(trainSet,trainClass,method,option);
    % predict
    [testClassPredicted,OtherOutput]=classificationPredict(model,testSet,testClass);
    % calculate the performance
    [classPerform,conMat]=perform(testClassPredicted,testClass,numel(unique([trainClass;testClass])));
    %     [testClassPredicted,classPerform,conMat,OtherOutput]=classification(trainSet,trainClass,testSet,testClass,method,option);
    testClassPredicteds(:,m)=testClassPredicted;
    classPerforms(m,:)=classPerform;
    conMats(:,:,m)=conMat;
    tElapseds(m,1)=OtherOutputTr{1}+OtherOutput{1};
    OtherOutputs{m}=OtherOutput;
end
end