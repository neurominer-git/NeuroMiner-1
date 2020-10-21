function [testClassPredicted,OtherOutput]=classificationPredict(model,testSet,testClass)
% predict the test samples by the model learned by "classificationTrain"
% Usage:
% [testClassPredicted,OtherOutput]=classificationPredict(model,testSet,testClass)
% testSet: matrix, the test set.
% testClass: column vector of numbers or string, the class labels of the
%       test/unknown set. It is actually unused in this function, thus, set it [].
% testClassPredicted: column vector, the predicted class labels of the test/unknown samples.
% OtherOutput: cell or row vector of cell.
%     If method=
%     'ksrsc': OtherOutput{1}=tElapsed; OtherOutput{2}=sparsity;
%     'subdic': OtherOutput{1}=tElapsed;
%     'svmm': OtherOutput{1}=tElapsed;
%     'hdlm': OtherOutput{1}=tElapsed;
%     'nmf', 'rnmf': OtherOutput{1}=tElapsed;
%     'svdd': OtherOutput{1}=tElapsed;
%     'hkm': OtherOutput{1}=tElapsed; OtherOutput{2}=beta;
%     'nc': OtherOutput{1}=tElapsed;
%     'lrc': OtherOutput{1}=tElapsed;
%     'logistic': OtherOutput{1}=tElapsed;
%     'bayesian': OtherOutput{1}=tElapsed; OtherOutput{2}=trainError; OtherOutput{3}=testPosterior;
%     'naiveBayes': OtherOutput{1}=tElapsed;
%     'knn': OtherOutput{1}=tElapsed;
%     'lsvm': OtherOutput{1}=tElapsed;

%%%%
% Copyright (C) <2012>  <Yifeng Li>
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
% 
% Contact Information:
% Yifeng Li
% University of Windsor
% li11112c@uwindsor.ca; yifeng.li.cn@gmail.com
% May 18, 2011
%%%%

option=model.trainOpt;

% downsample face images
if isfield(option,'ifDownSample')&& option.ifDownSample
    testSet=downsample(testSet,option.p);
end
% if tensor data
if option.ifTensor % only 3-order
    testSet=matrizicing(testSet,3);
    testSet=testSet';
end

% handle missing value, imputation
if option.ifMissValueImpute
    tfTest=isnan(testSet);
    ifMissValTestSet=any(any(tfTest));
    if ifMissValTestSet
        switch option.imputeMethod
            case 'knn'
                numTrain=size(model.trainSet,2);
                testSet=knnimpute([model.trainSet,testSet]);
                testSet=testSet(:,numTrain+1:end);
            case 'zero'
                testSet(tfTest)=0;
        end
    end
end

% normalization
if logical(option.normalization)
    switch option.normMethod
        case 'mean0std1'
            testSet=normmean0std1(testSet',option.trainSetMean,option.trainSetSTD);
            testSet=testSet';
        case 'unitl2norm'
            testSet=normc(testSet);
    end
end

% classification
switch model.method
    case 'naiveBayes'
        tic;
        testClassPredicted=model.modelBayes.predict(testSet');
        tElapsed=toc;
        OtherOutput{1}=tElapsed;
    case 'knn'
        tic;
        testClassPredicted=knnclassify(testSet',model.trainSet',model.trainClass,option.k);%KNN Classification
        tElapsed=toc;
        OtherOutput{1}=tElapsed;
    case 'nc'
        tic;
        [testClassPredicted,vals]=nearestCentroidPredict(model,testSet);
        tElapsed=toc;
%         perSVM=perform(testClassPredicted,testClass,model.numCl);
        OtherOutput{1}=tElapsed;        
%     case 'svm' % libsvm matlab toolbox is required
%         testSetting='-b 1';
%         testClass=zeros(size(testSet,2),1); % useless, just to fit the parameters
%         tic;
%         [predicted_label, accuracy, prob_estimates] = svmpredict(testClass, testSet', model.modelSVM,testSetting);
%         testClassPredicted=predicted_label;
%         tElapsed=toc;
%         OtherOutput{1}=tElapsed;
%         OtherOutput{2}= prob_estimates;
    case 'lssvm'
        tic;
        testClassPredicted=lssvmClassifier(trainSet,testSet,trainClass,option);
        tElapsed=toc;
        OtherOutput{1}=tElapsed;
%     case 'svmLi'
%         tic;
%         [testClassPredicted,vals]=softSVMPredict2(model,testSet);
%         tElapsed=toc;
% %         perSVM=perform(testClassPredicted,testClass,model.numCl);
%         OtherOutput{1}=tElapsed;
    case 'svmm'
        tic;
        [testClassPredicted,vals]=SVMMultiPredict(model,testSet);
        tElapsed=toc;
%         perSVM=perform(testClassPredicted,testClass,model.numCl);
        OtherOutput{1}=tElapsed;
    case 'svdd'
        tic;
        [testClassPredicted,vals]=SVDDPredict(model,testSet);
        tElapsed=toc;
%         perSVM=perform(testClassPredicted,testClass,model.numCl);
        OtherOutput{1}=tElapsed;
    case 'svddm'
        tic;
        [testClassPredicted,vals]=SVDDMPredict(model,testSet);
        tElapsed=toc;
%         perSVM=perform(testClassPredicted,testClass,model.numCl);
        OtherOutput{1}=tElapsed;
   case 'svdda'
        tic;
        [testClassPredicted,vals]=SVDDPredict(model,testSet);
        tElapsed=toc;
%         perSVM=perform(testClassPredicted,testClass,model.numCl);
        OtherOutput{1}=tElapsed;
    case 'lsvm' % local svm
        tic;
        testClassPredicted=localSVM(model.trainSet,model.trainClass,testSet,testClass,option);
        tElapsed=toc;
        OtherOutput{1}=tElapsed;
%     case 'kernelfda'
%         %kernel FDA
%         optionDefault.kernel='rbf';
%         optionDefault.kernelParameterValue=1;%1/size(trainSet,1);
%         option=mergeOption(option,optionDefault);
%         tic;
%         testClassPredicted=kfda(trainSet, testSet, trainClass, option.kernel,option.kernelParameterValue);
%         tElapsed=toc;
%         OtherOutput{1}=tElapsed;
%     case 'psdkernelfda'
%         optionDefault.numKernel=10;
%         optionDefault.lambda=10.^(-8);
%         log10Sigma=rand(optionDefault.numKernel,1)*3-1;
%         optionDefault.sigma=10.^(log10Sigma);
%         option=mergeOption(option,optionDefault);
%         tic;
%         [testClassPredicted,testValues,Theta]=OKernelFDA(trainSet,testSet,trainClass,option.lambda,option.sigma,option.numKernel,[]);
%         tElapsed=toc;
%         OtherOutput{1}=tElapsed;
%         OtherOutput{2}=Theta;
    case 'hkm'
        tic;
        testClassPredicted=hierarchialModelPredict(model,testSet);
        tElapsed=toc;
        OtherOutput{1}=tElapsed;
    case 'bayesian'
        optionDefault.type='diaglinear'; % 'linear','diaglinear','quadratic','diagquadratic','mahalanobis'
        option=mergeOption(option,optionDefault);
        tic;
        [testClassPredicted,trainError,testPosterior]=classify(testSet',trainSet',trainClass,option.type);
        tElapsed=toc;
        OtherOutput{1}=tElapsed;
        OtherOutput{2}=trainError;
        OtherOutput{3}=testPosterior;
    case 'nnls'
        tic;
        [testClassPredicted,sparsity]=nnlsClassifier(model.trainSet,model.trainClass,testSet,[],option);
        tElapsed=toc;
        OtherOutput{1}=tElapsed;
        OtherOutput{2}=sparsity;
    case 'bootstrapnnls'
        %        if size(trainSet,3)>1&& option.normalization
        %             cent=[0 0 1];
        %             scale=[0 0 0];
        %             [trainSet,Means,Scales]=nprocess(trainSet,cent,scale);
        %             testSet=nprocess(testSet,cent,scale,Means,Scales,1);
        %         end
        tic;
        testClassPredicted=bootstrapnnlsClassifier(model.trainSet,model.trainClass,testSet,testClass,option);
        tElapsed=toc;
        OtherOutput{1}=tElapsed;
%     case 'nnlslss'
%         tic;
%         if strcmp(option.kernel,'ds')
%             rank=1;
%             lambda=5;
%             option.param=[numR;numC;rank;lambda];
%         end
%         % model selection
%         if strcmp(option.kernel,'rbf') && option.search
%             option.normalization=0;
%             [param,maxAcc]=lineSearchNNLS(trainSet,trainClass,option);
%             option.param=param;
%         end
%         [testClassPredicted,sparsity]=nnlsClassifierLargeSampleSize(trainSet,trainClass,testSet,testClass,option);
%         tElapsed=toc;
%         OtherOutput{1}=tElapsed;
%         OtherOutput{2}=sparsity;
    case 'nmf'
        tic;
        [testClassPredicted]=nmfClassifier(model.trainSet,model.trainClass,testSet,testClass,option);
        tElapsed=toc;
        OtherOutput{1}=tElapsed;
    case 'rnmf'
        tic;
        [testClassPredicted]=repetitivenmfClassifier(model.trainSet,model.trainClass,testSet,testClass,option);
        tElapsed=toc;
        OtherOutput{1}=tElapsed;
    case 'src'
        tic;
        [testClassPredicted,sparsity]=src(trainSet,trainClass,testSet,option);
        tElapsed=toc;
        OtherOutput{1}=tElapsed;
        OtherOutput{2}=sparsity;
    case 'src2'
        tic;
        [testClassPredicted,sparsity]=SRC2(model.trainSet,model.trainClass,testSet,option);
        tElapsed=toc;
        OtherOutput{1}=tElapsed;
        OtherOutput{2}=sparsity;
    case 'ksrsc'
        tic;
        [testClassPredicted,sparsity]=KSRSCClassifier(model.trainSet,model.trainClass,testSet,option);
        tElapsed=toc;
        OtherOutput{1}=tElapsed;
        OtherOutput{2}=sparsity;
%     case 'ksrcl1ls'
%         tic;
%         optionDefault.ifModelSelection=false;
%         option=mergeOption(option,optionDefault);
%         if option.ifModelSelection
%             %             option.normalization=0;
%             %             [param,maxAcc]=lineSearchSRC2(trainSet,trainClass,option);
%             %             option.lambda=param;
%         end
%         [testClassPredicted, sparsity] = SRCl1LSKernel(trainSet, trainClass, testSet, option);
%         tElapsed=toc;
%         OtherOutput{1}=tElapsed;
%         OtherOutput{2}=sparsity;
%     case 'ksrcnnls'
%         tic;
%         optionDefault.ifModelSelection=false;
%         option=mergeOption(option,optionDefault);
%         if option.ifModelSelection
%             %             option.normalization=0;
%             %             [param,maxAcc]=lineSearchSRC2(trainSet,trainClass,option);
%             %             option.lambda=param;
%         end
%         [testClassPredicted, sparsity] = SRCNNLSKernel(trainSet, trainClass, testSet, option);
%         tElapsed=toc;
%         OtherOutput{1}=tElapsed;
%         OtherOutput{2}=sparsity;
%     case 'ksrcl1nnls'
%         tic;
%         optionDefault.ifModelSelection=false;
%         option=mergeOption(option,optionDefault);
%         if option.ifModelSelection
%             %             option.normalization=0;
%             %             [param,maxAcc]=lineSearchSRC2(trainSet,trainClass,option);
%             %             option.lambda=param;
%         end
%         [testClassPredicted, sparsity] = SRCl1NNLSKernel(trainSet, trainClass, testSet, option);
%         tElapsed=toc;
%         OtherOutput{1}=tElapsed;
%         OtherOutput{2}=sparsity;
%     case 'ksrcl1lsdl'
%         tic;
%         optionDefault.ifModelSelection=false;
%         optionDefault.dictLearn=true;
%         optionDefault.k=10;
%         optionDefault.classifierAfterSR='hdlm';
%         optionDefault.classifierAfterSRKernel='rbf';
%         optionDefault.classifierAfterSRParam=2^0;
%         option=mergeOption(option,optionDefault);
%         option.SRMethod='l1ls';
%         if option.ifModelSelection
%             %             option.normalization=0;
%             %             [param,maxAcc]=lineSearchSRC2(trainSet,trainClass,option);
%             %             option.lambda=param;
%         end
%         
%         if option.dictLearn
%             % dictionary learning
%             trainTrain=computeKernelMatrix(trainSet,trainSet,option);
%             trainTest=computeKernelMatrix(trainSet,testSet,option);
%             testTest=computeKernelMatrix(testSet,testSet,option);
%             testTest=diag(testTest);
%             [AtA,Y,~,~,~]=KSRDL(trainTrain,option.k,option);
%             AtTest=Y'\trainTest;
%             % sparse coding
%             [X,~,sparsity]=KSRSC(AtA,AtTest,testTest,option);
%             % use a classifier
%             optCl=option;
%             [testClassPredicted,~,~]=classification(Y,trainClass,X,testClass,option.classifierAfterSR,optCl);
%         else % do not learn dictionary
%             option.classifierAfterSR='ksrcl1ls';
% %             option.classifierAfterSRKernel='rbf';
% %             option.classifierAfterSRParam=2^0;
%             [testClassPredicted, sparsity] = SRCl1LSKernel(trainSet, trainClass, testSet, option);
%         end
%         tElapsed=toc;
%         OtherOutput{1}=tElapsed;
%         OtherOutput{2}=sparsity;
%     case 'ksrcnnlsdl'
%         tic;
%         optionDefault.ifModelSelection=false;
%         optionDefault.dictLearn=true;
%         optionDefault.k=10;
%         optionDefault.classifierAfterSR='hdlm';
%         optionDefault.classifierAfterSRKernel='rbf';
%         optionDefault.classifierAfterSRParam=2^0;
%         option=mergeOption(option,optionDefault);
%         option.SRMethod='nnls';
%         if option.ifModelSelection
%             %             option.normalization=0;
%             %             [param,maxAcc]=lineSearchSRC2(trainSet,trainClass,option);
%             %             option.lambda=param;
%         end
%         
%         if option.dictLearn
%             % dictionary learning
%             trainTrain=computeKernelMatrix(trainSet,trainSet,option);
%             trainTest=computeKernelMatrix(trainSet,testSet,option);
%             testTest=computeKernelMatrix(testSet,testSet,option);
%             testTest=diag(testTest);
%             [AtA,Y,~,~,~]=KSRDL(trainTrain,option.k,option);
%             AtTest=Y'\trainTest;
%             % sparse coding
%             [X,~,sparsity]=KSRSC(AtA,AtTest,testTest,option);
%             % use a classifier
%             optCl=option;
%             [testClassPredicted,~,~]=classification(Y,trainClass,X,testClass,option.classifierAfterSR,optCl);
%         else % do not learn dictionary
%             option.classifierAfterSR='ksrcnnls';
% %             option.classifierAfterSRKernel='rbf';
% %             option.classifierAfterSRParam=2^0;
%             [testClassPredicted, sparsity] = SRCl1LSKernel(trainSet, trainClass, testSet, option);
%         end
%         tElapsed=toc;
%         OtherOutput{1}=tElapsed;
%         OtherOutput{2}=sparsity;
%     case 'ksrcl1nnlsdl'
%         tic;
%         optionDefault.ifModelSelection=false;
%         optionDefault.dictLearn=true;
%         optionDefault.k=10;
%         optionDefault.classifierAfterSR='hdlm';
%         optionDefault.classifierAfterSRKernel='rbf';
%         optionDefault.classifierAfterSRParam=2^0;
%         option=mergeOption(option,optionDefault);
%         option.SRMethod='l1nnls';
%         if option.ifModelSelection
%             %             option.normalization=0;
%             %             [param,maxAcc]=lineSearchSRC2(trainSet,trainClass,option);
%             %             option.lambda=param;
%         end
%         
%         if option.dictLearn
%             % dictionary learning
%             trainTrain=computeKernelMatrix(trainSet,trainSet,option);
%             trainTest=computeKernelMatrix(trainSet,testSet,option);
%             testTest=computeKernelMatrix(testSet,testSet,option);
%             testTest=diag(testTest);
%             [AtA,Y,~,~,~]=KSRDL(trainTrain,option.k,option);
%             AtTest=Y'\trainTest;
%             % sparse coding
%             [X,~,sparsity]=KSRSC(AtA,AtTest,testTest,option);
%             % use a classifier
%             optCl=option;
%             [testClassPredicted,~,~]=classification(Y,trainClass,X,testClass,option.classifierAfterSR,optCl);
%         else % do not learn dictionary
%             option.classifierAfterSR='ksrcl1nnls';
% %             option.classifierAfterSRKernel='rbf';
% %             option.classifierAfterSRParam=2^0;
%             [testClassPredicted, sparsity] = SRCl1LSKernel(trainSet, trainClass, testSet, option);
%         end
%         tElapsed=toc;
%         OtherOutput{1}=tElapsed;
%         OtherOutput{2}=sparsity;
    case 'msrc'
        tic;
        [testClassPredicted, sparsity] = SRC2(model.metaSample, model.metaClass, testSet, option);
        tElapsed=toc;
        OtherOutput{1}=tElapsed;
        OtherOutput{2}=sparsity;
    case 'mnnls'
        tic;
        optionDefault.ifModelSelection=false;
        optionDefault.metaSampleMethod='nmf';
        option=mergeOption(option,optionDefault);
        [trainSet,trainClass]=computeMetaSample(trainSet,trainClass,option);
        [testClassPredicted,sparsity]=nnlsClassifier(trainSet,trainClass,testSet,testClass,option);
%         methodCl='svm'; % SVM does not work
%         optionsvm.normalization=0;
%         optionsvm.trainSetting='-t 2 -c 1 -b 1';
%         optionsvm.testSetting='-b 1';
%         optionsvm.search=false;
%         [testClassPredicted,classPerform,OtherOutput]=classification(trainSet,trainClass,testSet,testClass,methodCl,optionsvm);
%         sparsity=[];
        tElapsed=toc;
        OtherOutput{1}=tElapsed;
        OtherOutput{2}=sparsity;
    case 'lrc'
        tic;
        [testClassPredicted]=lrc(model.trainSet,model.trainClass,testSet,testClass,option);
        tElapsed=toc;
        OtherOutput{1}=tElapsed;
%     case 'ls'
%         tic;
%         option.searchThreshold=false;
%         [testClassPredicted,beta]=lsClassifier(trainSet,testSet,trainClass,testClass,option);
%         %         marginLS=2/sqrt(sum(beta.^2))
%         %         trainLabels=(beta'*[ones(1,size(trainSet,2));trainSet])';
%         tElapsed=toc;
%         OtherOutput{1}=tElapsed;
%         OtherOutput{2}=beta;
%     case 'lsm'
%         tic;
%         [testClassPredicted,beta]=lsmClassifier(trainSet,testSet,trainClass,testClass);
%         tElapsed=toc;
%         OtherOutput{1}=tElapsed;
%         OtherOutput{2}=beta;
%     case 'lsmBag'
%         optionDefault.numRandom=99;
%         option=mergeOption(option,optionDefault);
%         tic;
%         testClassPredicted=bootstrapls(trainSet,trainClass,testSet,testClass,option);
%         tElapsed=toc;
%         OtherOutput{1}=tElapsed;
    case 'hdlm'
         tic;
        testClassPredicted=khdlmPredict(model,testSet);
        tElapsed=toc;
%         perSVM=perform(testClassPredicted,testClass,model.numCl);
        OtherOutput{1}=tElapsed;
        tic;
    case 'elm'
        tic;
        Elm_Type=1;
        [testClassPredicted, TrainingTime, TestingTime, TrainingAccuracy, TestingAccuracy] = ...
            elmModified(model.trainSet, testSet, model.trainClass, testClass ,Elm_Type, option.NumberofHiddenNeurons, option.ActivationFunction);
        tElapsed=toc;
        OtherOutput{1}=tElapsed;
%     case 'elmBag'
%         tic;
%         optionDefault.NumberofHiddenNeurons=200;
%         optionDefault.ActivationFunction='sig';
%         option=mergeOption(option,optionDefault);
%         testClassPredicted = bootstrapelm(trainSet,trainClass,testSet,testClass, option);
%         tElapsed=toc;
%         OtherOutput{1}=tElapsed;
    case 'logistic'
        tic;
        testClassPredicted=logisticRegressPredict(model,testSet);
        tElapsed=toc;
        OtherOutput{1}=tElapsed;
%     case 'naive'
%         tic;
%         nb = NaiveBayes.fit(trainSet', trainClass);
%         testClassPredicted = predict(nb,testSet');
%         tElapsed=toc;
%         OtherOutput{1}=tElapsed;
%     case 'tree'
%         tic;
%         t = classregtree(trainSet',trainClass,'Method','classification');
%         testClassPredicted  = eval(t,testSet');
%         testClassPredicted= cell2mat(testClassPredicted);
%         tElapsed=toc;
%         OtherOutput{1}=tElapsed;
    otherwise
        error('Please a correct method parameter for a classifier.');
end
end