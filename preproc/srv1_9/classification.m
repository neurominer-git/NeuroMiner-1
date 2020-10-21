function [testClassPredicted,classPerform,conMat,OtherOutput]=classification(trainSet,trainClass,testSet,testClass,method,option)
% Classification using NNLS, BNNLS, NMF, RNMF, KNN, SVM, Bayesian, Kernel-FDA (for future extension), or Optimal-Kernel-FDA (for future extension) classifiers.
% Usage:
% [testClassPredicted,classPerform,OtherOutput]=classification(trainSet,trainClass,[],testClass)
% [testClassPredicted,sparsity]=classification(trainSet,trainClass,testSet,testClass)
% [testClassPredicted,sparsity]=classification(trainSet,trainClass,testSet,testClass,option)
% trainSet, matrix, the training set with samples in columns and features in rows.
% trainClass: column vector of numbers or string, the class labels of the traning set.
% testSet: matrix, the test set.
% testClass: column vector of numbers or string, the class labels of the
%       test/unknown set. It is actually unused in this function, thus, set it [].
% method: string, designate which classifier to use. It could be
%     'nnls': the NNLS classifier;
%     'bootstrapnnls': the Bootstrap based NNLS classifier;
%     'nmf': the NMF classifier;
%     'rnmf': the RNMF classifier;
%     'knn': the KNN classifier;
%     'svm': the SVM classifier;
%     'bayesian': the Bayesian classifier.
%%%     'kernelfda': (for future extension) the Kernel Fisher Discriminative Analysis classifier;
%%%    'psdkernelfda': (for future extension) the Optimal Kernel Fisher Discriminative Analysis classifier;
% option: struct: the options to configue specific classifier:
%     If method=
%     'nnls':
%            option.normalization, scalar, 0: no normaization (default), 1: normalize the data under each feature to have mean 1 and std 1;
%            The rest fields of option is the same as the input option for "nnlsClassifier" function. Type "help nnlsClassifier" for more information.
%     'bootstrapnnls':
%            option.normalization, scalar, 0: no normaization (default), 1: normalize the data under each feature to have mean 1 and std 1;
%            The rest fields of option is the same as the input option for "bootstrapnnlsClassifier" function. Type "help bootstrapnnlsClassifier" for more information.
%     'nmf':
%            option.normalization, scalar, 0: no normaization (default), 1: normalize the data under each feature to have mean 1 and std 1;
%            The rest fields of option is the same as the input option for "nmfClassifier" function. Type "help nmfClassifier" for more information.
%     'rnmf':
%            option.normalization, scalar, 0: no normaization (default), 1: normalize the data under each feature to have mean 1 and std 1;
%            The rest fields of option is the same as the input option for "repetitivenmfClassifier" function. Type "help repetitivenmfClassifier" for more information.
%     'knn':
%            option.normalization: scalar, 0: no normaization (default), 1: normalize the data under each feature to have mean 1 and std 1;
%            option.k: the number of neares neighbors. The default is 1.
%     'svm':
%            option.normalization: scalar, 0: no normaization (default), 1: normalize the data under each feature to have mean 1 and std 1;
%            option.trainSetting: string, the options for the "svmtrain" function. The default is '-s 0 -t 2 -b 1'. See the LIBSVM: http://www.csie.ntu.edu.tw/~cjlin/libsvm;
%            option.testSetting: string, the options for the "svmpredict" function. The default is '-b 1'. See the LIBSVM: http://www.csie.ntu.edu.tw/~cjlin/libsvm.
%     'bayesian':
%            option.normalization: scalar, 0: no normaization (default), 1: normalize the data under each feature to have mean 1 and std 1;
%            option.type: string, specify the type of discriminant function.
%                 It could be 'linear','diaglinear','quadratic','diagquadratic','mahalanobis'. Type "help classify" for more information.
%%%     'kernelfda': (for future extension)
%            option.normalization: scalar, 0: no normaization (default), 1: normalize the data under each feature to have mean 1 and std 1;
%            option.kernel: string, the kernel function. Type "kfda" for more information. It could be
%                  'linear':  option.param=[];
%                  'polynomial': option.param is [Gamma;Coefficient;Degree], the default is [1;0;2];
%                  'rbf': option.param is sigma, the default is 1/#features;
%                  'sigmoid': option.param is [alpha;beta], the default is [1;0];
%                  you own kernel function name.
%%%     'psdkernelfda': (for future extension) the Optimal Kernel Fisher Discriminative Analysis classifier;
% testClassPredicted: column vector, the predicted class labels of the test/unknown samples.
% classPerform: row vector, the classification performance.
%     If binary-class problem, classPerform includes PPV, NPV, Specificity, Sensitivity, Accuracy, and Balanced Accuracy.
%     If multi-class problem, classPerform includes accuracy for class 1, accuracy of class 2, ..., accuracy, balanced accuracy.
% OtherOutput: cell or row vector of cell.
%     If method=
%     'nnls': OtherOutput{1}=tElapsed; OtherOutput{2}=sparsity;
%     'bootstrapnnls': OtherOutput{1}=tElapsed;
%     'nmf': OtherOutput{1}=tElapsed;
%     'rnmf': OtherOutput{1}=tElapsed;
%     'knn': OtherOutput{1}=tElapsed;
%     'svm': OtherOutput{1}=tElapsed; OtherOutput{2}= prob_estimates;
%     'bayesian': OtherOutput{1}=tElapsed; OtherOutput{2}=trainError; OtherOutput{3}=testPosterior;
% Contact Information:
% Yifeng Li
% University of Windsor
% li11112c@uwindsor.ca; yifeng.li.cn@gmail.com
% May 18, 2011

if nargin<6
    option=[];
end
OtherOutput=[];
optionDefault.normalization=0;
optionDefault.normMethod='mean0std1';
optionDefault.ifDownSample=false;
optionDefault.p=16; % downsample rate
optionDefault.ifMissValueImpute=true;
optionDefault.imputeMethod='knn';
optionDefault.search=false;
option=mergeOption(option,optionDefault);

% downsample face images
if isfield(option,'ifDownSample')&& option.ifDownSample
    trainSet=downsample(trainSet,option.p);
    testSet=downsample(testSet,option.p);
end
% if tensor data
if size(trainSet,3)>1
    numR=size(trainSet,1);
    numC=size(trainSet,2);
    trainSet=matrizicing(trainSet,3);
    testSet=matrizicing(testSet,3);
    trainSet=trainSet';
    testSet=testSet';
end

% handle missing value, imputation
if option.ifMissValueImpute
    tfTrain=isnan(trainSet);
    tfTest=isnan(testSet);
    ifMissValTrainSet=any(any(tfTrain));
    ifMissValTestSet=any(any(tfTest));
    if ifMissValTrainSet
        switch option.imputeMethod
            case 'knn'
                trainSet=knnimpute(trainSet);
            case 'zero'
                trainSet(tfTrain)=0;
        end
    end
    if ifMissValTestSet
        switch option.imputeMethod
            case 'knn'
                numTrain=size(trainSet,2);
                testSet=knnimpute([trainSet,testSet]);
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
            [trainSet,trainSetMean,trainSetSTD]=normmean0std1(trainSet');
            trainSet=trainSet';
            testSet=normmean0std1(testSet',trainSetMean,trainSetSTD);
            testSet=testSet';
        case 'unitl2norm'
            trainSet=normc(trainSet);
            testSet=normc(testSet);
    end
end

% TFTr=isnan(trainSet);
% if any(TFTr)
%     trainSet(TFTr)=0;
% end
% TFTe=isnan(testSet);
% if any(TFTe)
%     testSet(TFTr)=0;
% end

% classification
switch method
    case 'knn'
        optionDefault.k=1;
        option=mergeOption(option,optionDefault);
        tic;
        testClassPredicted=knnclassify(testSet',trainSet',trainClass,option.k);%KNN Classification
        tElapsed=toc;
        OtherOutput{1}=tElapsed;
    case 'svm' % libsvm matlab toolbox is required
        optionDefault.trainSetting='-s 0 -t 2 -b 1';
        optionDefault.testSetting='-b 1';
        optionDefault.search=false;
        option=mergeOption(option,optionDefault);
        tic;
        % model selection
        if (isempty(strfind(option.trainSetting,'-t'))|| ~isempty(strfind(option.trainSetting,'-t 2'))) && option.search % rbf and search
            option.normalization=0;
            [gammaOptimal,COptimal,accMax]=gridSearch(trainSet,trainClass,option);
            option.trainSetting=['-s 0 -t 2 -g ',num2str(gammaOptimal),' -c ',num2str(COptimal),'  -b 1'];
        end
        if (isempty(strfind(option.trainSetting,'-t'))|| ~isempty(strfind(option.trainSetting,'-t 0'))) && option.search % linear and search
            option.normalization=0;
            [COptimal,accMax]=linearSearchSVM(trainSet,trainClass,option);
            option.trainSetting=['-s 0 -t 0 -c ',num2str(COptimal),'  -b 1'];
        end
        model=svmtrain(trainClass,trainSet',option.trainSetting);
        % margin, use hdlm to compute
        %         betaSVM=pinv([ones(model.totalSV,1),full(model.SVs)])*[-ones(model.nSV(1),1);ones(model.nSV(2),1)];
        %         marginSVM=2/sqrt(sum(betaSVM.^2));
        marginSVM=0;
        %         trainLabelsSVM=(betaSVM'*[ones(1,size(trainSet,2));trainSet])';
        [predicted_label, accuracy, prob_estimates] = svmpredict(testClass, testSet', model,option.testSetting);
        testClassPredicted=predicted_label;
        tElapsed=toc;
        OtherOutput{1}=tElapsed;
        OtherOutput{2}= prob_estimates;
        OtherOutput{3}=marginSVM;
    case 'lssvm'
        tic;
        testClassPredicted=lssvmClassifier(trainSet,testSet,trainClass,option);
        tElapsed=toc;
        OtherOutput{1}=tElapsed;
    case 'svmLi'
        tic;
        % model selection
        if option.search
            optGS.classifier='svmLi';
            [alphaOptimal,betaOptimal,accMax]=gridSearchUniverse(trainSet,trainClass,optGS);
            option.C=alphaOptimal;
            option.param=betaOptimal;
        end
        model=softSVMTrain2(trainSet,trainClass,option);
        [testClassPredicted,vals]=softSVMPredict2(model,testSet);
        tElapsed=toc;
%         perSVM=perform(testClassPredicted,testClass,model.numCl);
        OtherOutput{1}=tElapsed;
    case 'lsvm' % local svm
        tic;
        testClassPredicted=localSVM(trainSet,trainClass,testSet,testClass,option);
        tElapsed=toc;
        OtherOutput{1}=tElapsed;
    case 'kernelfda'
        %kernel FDA
        optionDefault.kernel='rbf';
        optionDefault.kernelParameterValue=1;%1/size(trainSet,1);
        option=mergeOption(option,optionDefault);
        tic;
        testClassPredicted=kfda(trainSet, testSet, trainClass, option.kernel,option.kernelParameterValue);
        tElapsed=toc;
        OtherOutput{1}=tElapsed;
    case 'psdkernelfda'
        optionDefault.numKernel=10;
        optionDefault.lambda=10.^(-8);
        log10Sigma=rand(optionDefault.numKernel,1)*3-1;
        optionDefault.sigma=10.^(log10Sigma);
        option=mergeOption(option,optionDefault);
        tic;
        [testClassPredicted,testValues,Theta]=OKernelFDA(trainSet,testSet,trainClass,option.lambda,option.sigma,option.numKernel,[]);
        tElapsed=toc;
        OtherOutput{1}=tElapsed;
        OtherOutput{2}=Theta;
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
        if strcmp(option.kernel,'ds')
            rank=1;
            lambda=5;
            option.param=[numR;numC;rank;lambda];
        end
        % model selection
        if strcmp(option.kernel,'rbf') && option.search
            option.normalization=0;
            [param,maxAcc]=lineSearchNNLS(trainSet,trainClass,option);
            option.param=param;
        end
        [testClassPredicted,sparsity]=nnlsClassifier(trainSet,trainClass,testSet,testClass,option);
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
        optionDefault.method='nnls';%
        optionDefault.predicter='max';
        optionDefault.k=1;
        optionDefault.numRandom=99;
        option=mergeOption(option,optionDefault);
        tic;
        if strcmp(option.kernel,'ds')
            rank=NaN;
            lambda=NaN;
            option.param=[numR;numC;rank;lambda];
        end
        testClassPredicted=bootstrapnnlsClassifier(trainSet,trainClass,testSet,testClass,option);
        tElapsed=toc;
        OtherOutput{1}=tElapsed;
    case 'nnlslss'
        tic;
        if strcmp(option.kernel,'ds')
            rank=1;
            lambda=5;
            option.param=[numR;numC;rank;lambda];
        end
        % model selection
        if strcmp(option.kernel,'rbf') && option.search
            option.normalization=0;
            [param,maxAcc]=lineSearchNNLS(trainSet,trainClass,option);
            option.param=param;
        end
        [testClassPredicted,sparsity]=nnlsClassifierLargeSampleSize(trainSet,trainClass,testSet,testClass,option);
        tElapsed=toc;
        OtherOutput{1}=tElapsed;
        OtherOutput{2}=sparsity;
    case 'nmf'
        optionDefault.facts=numel(unique([trainClass;testClass]));
        option=mergeOption(option,optionDefault);
        tic;
        [testClassPredicted]=nmfClassifier(trainSet,trainClass,testSet,testClass,option);
        tElapsed=toc;
        OtherOutput{1}=tElapsed;
    case 'rnmf'
        optionDefault.facts=7;%numel(unique([trainClass;testClass]));
        option=mergeOption(option,optionDefault);
        tic;
        [testClassPredicted]=repetitivenmfClassifier(trainSet,trainClass,testSet,testClass,option);
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
        optionDefault.ifModelSelection=false;
        option=mergeOption(option,optionDefault);
        if option.ifModelSelection
            option.normalization=0;
            [param,maxAcc]=lineSearchSRC2(trainSet,trainClass,option);
            option.lambda=param;
        end
        [testClassPredicted, sparsity] = SRC2(trainSet, trainClass, testSet, option);
        tElapsed=toc;
        OtherOutput{1}=tElapsed;
        OtherOutput{2}=sparsity;
    case 'ksrcl1ls'
        tic;
        optionDefault.ifModelSelection=false;
        option=mergeOption(option,optionDefault);
        if option.ifModelSelection
            %             option.normalization=0;
            %             [param,maxAcc]=lineSearchSRC2(trainSet,trainClass,option);
            %             option.lambda=param;
        end
        [testClassPredicted, sparsity] = SRCl1LSKernel(trainSet, trainClass, testSet, option);
        tElapsed=toc;
        OtherOutput{1}=tElapsed;
        OtherOutput{2}=sparsity;
    case 'ksrcnnls'
        tic;
        optionDefault.ifModelSelection=false;
        option=mergeOption(option,optionDefault);
        if option.ifModelSelection
            %             option.normalization=0;
            %             [param,maxAcc]=lineSearchSRC2(trainSet,trainClass,option);
            %             option.lambda=param;
        end
        [testClassPredicted, sparsity] = SRCNNLSKernel(trainSet, trainClass, testSet, option);
        tElapsed=toc;
        OtherOutput{1}=tElapsed;
        OtherOutput{2}=sparsity;
    case 'ksrcl1nnls'
        tic;
        optionDefault.ifModelSelection=false;
        option=mergeOption(option,optionDefault);
        if option.ifModelSelection
            %             option.normalization=0;
            %             [param,maxAcc]=lineSearchSRC2(trainSet,trainClass,option);
            %             option.lambda=param;
        end
        [testClassPredicted, sparsity] = SRCl1NNLSKernel(trainSet, trainClass, testSet, option);
        tElapsed=toc;
        OtherOutput{1}=tElapsed;
        OtherOutput{2}=sparsity;
    case 'ksrcl1lsdl'
        tic;
        optionDefault.ifModelSelection=false;
        optionDefault.dictLearn=true;
        optionDefault.k=10;
        optionDefault.classifierAfterSR='hdlm';
        optionDefault.classifierAfterSRKernel='rbf';
        optionDefault.classifierAfterSRParam=2^0;
        option=mergeOption(option,optionDefault);
        option.SRMethod='l1ls';
        if option.ifModelSelection
            %             option.normalization=0;
            %             [param,maxAcc]=lineSearchSRC2(trainSet,trainClass,option);
            %             option.lambda=param;
        end
        
        if option.dictLearn
            % dictionary learning
            trainTrain=computeKernelMatrix(trainSet,trainSet,option);
            trainTest=computeKernelMatrix(trainSet,testSet,option);
            testTest=computeKernelMatrix(testSet,testSet,option);
            testTest=diag(testTest);
            [AtA,Y,~,~,~]=KSRDL(trainTrain,option.k,option);
            AtTest=Y'\trainTest;
            % sparse coding
            [X,~,sparsity]=KSRSC(AtA,AtTest,testTest,option);
            % use a classifier
            optCl=option;
            [testClassPredicted,~,~]=classification(Y,trainClass,X,testClass,option.classifierAfterSR,optCl);
        else % do not learn dictionary
            option.classifierAfterSR='ksrcl1ls';
%             option.classifierAfterSRKernel='rbf';
%             option.classifierAfterSRParam=2^0;
            [testClassPredicted, sparsity] = SRCl1LSKernel(trainSet, trainClass, testSet, option);
        end
        tElapsed=toc;
        OtherOutput{1}=tElapsed;
        OtherOutput{2}=sparsity;
    case 'ksrcnnlsdl'
        tic;
        optionDefault.ifModelSelection=false;
        optionDefault.dictLearn=true;
        optionDefault.k=10;
        optionDefault.classifierAfterSR='hdlm';
        optionDefault.classifierAfterSRKernel='rbf';
        optionDefault.classifierAfterSRParam=2^0;
        option=mergeOption(option,optionDefault);
        option.SRMethod='nnls';
        if option.ifModelSelection
            %             option.normalization=0;
            %             [param,maxAcc]=lineSearchSRC2(trainSet,trainClass,option);
            %             option.lambda=param;
        end
        
        if option.dictLearn
            % dictionary learning
            trainTrain=computeKernelMatrix(trainSet,trainSet,option);
            trainTest=computeKernelMatrix(trainSet,testSet,option);
            testTest=computeKernelMatrix(testSet,testSet,option);
            testTest=diag(testTest);
            [AtA,Y,~,~,~]=KSRDL(trainTrain,option.k,option);
            AtTest=Y'\trainTest;
            % sparse coding
            [X,~,sparsity]=KSRSC(AtA,AtTest,testTest,option);
            % use a classifier
            optCl=option;
            [testClassPredicted,~,~]=classification(Y,trainClass,X,testClass,option.classifierAfterSR,optCl);
        else % do not learn dictionary
            option.classifierAfterSR='ksrcnnls';
%             option.classifierAfterSRKernel='rbf';
%             option.classifierAfterSRParam=2^0;
            [testClassPredicted, sparsity] = SRCl1LSKernel(trainSet, trainClass, testSet, option);
        end
        tElapsed=toc;
        OtherOutput{1}=tElapsed;
        OtherOutput{2}=sparsity;
    case 'ksrcl1nnlsdl'
        tic;
        optionDefault.ifModelSelection=false;
        optionDefault.dictLearn=true;
        optionDefault.k=10;
        optionDefault.classifierAfterSR='hdlm';
        optionDefault.classifierAfterSRKernel='rbf';
        optionDefault.classifierAfterSRParam=2^0;
        option=mergeOption(option,optionDefault);
        option.SRMethod='l1nnls';
        if option.ifModelSelection
            %             option.normalization=0;
            %             [param,maxAcc]=lineSearchSRC2(trainSet,trainClass,option);
            %             option.lambda=param;
        end
        
        if option.dictLearn
            % dictionary learning
            trainTrain=computeKernelMatrix(trainSet,trainSet,option);
            trainTest=computeKernelMatrix(trainSet,testSet,option);
            testTest=computeKernelMatrix(testSet,testSet,option);
            testTest=diag(testTest);
            [AtA,Y,~,~,~]=KSRDL(trainTrain,option.k,option);
            AtTest=Y'\trainTest;
            % sparse coding
            [X,~,sparsity]=KSRSC(AtA,AtTest,testTest,option);
            % use a classifier
            optCl=option;
            [testClassPredicted,~,~]=classification(Y,trainClass,X,testClass,option.classifierAfterSR,optCl);
        else % do not learn dictionary
            option.classifierAfterSR='ksrcl1nnls';
%             option.classifierAfterSRKernel='rbf';
%             option.classifierAfterSRParam=2^0;
            [testClassPredicted, sparsity] = SRCl1LSKernel(trainSet, trainClass, testSet, option);
        end
        tElapsed=toc;
        OtherOutput{1}=tElapsed;
        OtherOutput{2}=sparsity;
    case 'msrc'
        tic;
        optionDefault.ifModelSelection=false;
        optionDefault.metaSampleMethod='svd';
        option=mergeOption(option,optionDefault);
        [trainSet,trainClass]=computeMetaSample(trainSet,trainClass,option);
        if option.ifModelSelection
            option.normalization=0;
            [param,maxAcc]=lineSearchSRC2(trainSet,trainClass,option);
            option.lambda=param;
        end
        [testClassPredicted, sparsity] = SRC2(trainSet, trainClass, testSet, option);
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
        [testClassPredicted]=lrc(trainSet,trainClass,testSet,testClass,option);
        tElapsed=toc;
        OtherOutput{1}=tElapsed;
    case 'ls'
        tic;
        option.searchThreshold=false;
        [testClassPredicted,beta]=lsClassifier(trainSet,testSet,trainClass,testClass,option);
        %         marginLS=2/sqrt(sum(beta.^2))
        %         trainLabels=(beta'*[ones(1,size(trainSet,2));trainSet])';
        tElapsed=toc;
        OtherOutput{1}=tElapsed;
        OtherOutput{2}=beta;
    case 'lsm'
        tic;
        [testClassPredicted,beta]=lsmClassifier(trainSet,testSet,trainClass,testClass);
        tElapsed=toc;
        OtherOutput{1}=tElapsed;
        OtherOutput{2}=beta;
    case 'lsmBag'
        optionDefault.numRandom=99;
        option=mergeOption(option,optionDefault);
        tic;
        testClassPredicted=bootstrapls(trainSet,trainClass,testSet,testClass,option);
        tElapsed=toc;
        OtherOutput{1}=tElapsed;
    case 'hdlm'
        tic;
        % model selection
        if strcmp(option.kernel,'rbf') && option.search
            option.normalization=0;
            [param,maxAcc]=lineSearchHDLM(trainSet,trainClass,option);
            option.param=param;
        end
        if strcmp(option.kernel,'sigmoid') && option.search
            option.normalization=0;
            [beta,alpha,maxAcc]=gridSearchHDLM(trainSet,trainClass,option);
            option.param=[alpha,beta];
        end
        [testClassPredicted,beta]=khdlmClassifier(trainSet,testSet,trainClass,testClass,option);
        % margin
        %         marginHDLM=2/sqrt(sum(beta.^2));
        marginHDLM=0;
        tElapsed=toc;
        OtherOutput{1}=tElapsed;
        OtherOutput{2}=beta;
        OtherOutput{3}=marginHDLM;
    case 'elm'
        tic;
        Elm_Type=1;
        optionDefault.TimesofHiddenNeurons=10;
        optionDefault.ActivationFunction='sig';
        option=mergeOption(option,optionDefault);
        option.NumberofHiddenNeurons=option.TimesofHiddenNeurons * size(trainSet,2);
        %         if option.NumberofHiddenNeurons==0
        %             option.NumberofHiddenNeurons=size(trainSet,2);
        %         end
        [testClassPredicted, TrainingTime, TestingTime, TrainingAccuracy, TestingAccuracy] = ...
            elmModified(trainSet, testSet, trainClass, testClass ,Elm_Type, option.NumberofHiddenNeurons, option.ActivationFunction);
        tElapsed=toc;
        OtherOutput{1}=tElapsed;
    case 'elmBag'
        tic;
        optionDefault.NumberofHiddenNeurons=200;
        optionDefault.ActivationFunction='sig';
        option=mergeOption(option,optionDefault);
        testClassPredicted = bootstrapelm(trainSet,trainClass,testSet,testClass, option);
        tElapsed=toc;
        OtherOutput{1}=tElapsed;
    case 'logistic'
        tic;
        optionDefault=[];
        option=mergeOption(option,optionDefault);
        model=logisticRegressTrain(trainSet,trainClass,option);
        testClassPredicted=logisticRegressPredict(model,testSet);
        tElapsed=toc;
        OtherOutput{1}=tElapsed;
    case 'naive'
        tic;
        nb = NaiveBayes.fit(trainSet', trainClass);
        testClassPredicted = predict(nb,testSet');
        tElapsed=toc;
        OtherOutput{1}=tElapsed;
    case 'tree'
        tic;
        t = classregtree(trainSet',trainClass,'Method','classification');
        testClassPredicted  = eval(t,testSet');
        testClassPredicted= cell2mat(testClassPredicted);
        tElapsed=toc;
        OtherOutput{1}=tElapsed;
    otherwise
        error('Please a correct method parameter for a classifier.');
end
% calculate the performance
[classPerform,conMat]=perform(testClassPredicted,testClass,numel(unique([trainClass;testClass])));
end