function [model,OtherOutput]=classificationTrain(trainSet,trainClass,method,option)
% Learn a classification model specified by "method".
% Usage:
% [model,OtherOutput]=classificationTrain(trainSet,trainClass,method)
% [model,OtherOutput]=classificationTrain(trainSet,trainClass,method,option)
% trainSet, matrix, the training set with samples in columns and features in rows.
% trainClass: column vector of numbers or string, the class labels of the traning set.
% method: string, designate which classifier to use. It could be
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
% option: struct: the options to configue specific classifier:
% option.normalization: scalar, 0: no normaization (default), 1: normalize the data under each feature to have mean 1 and std 1;
% For the rest fields of option, please see the corresponding classifier.
% if method =='knn', option.k: the number of neares neighbors. The default is 1.
% if method=='bayesian',  option.type: string, specify the type of discriminant function. It could be 'linear','diaglinear','quadratic','diagquadratic','mahalanobis'. Type "help classify" for more information.
% model: struct, the specific model learned.
% OtherOutput: cell or row vector of cell.
%     If method=
%     'ksrsc': OtherOutput{1}=tElapsed; % time elapsed.
%     'subdic': OtherOutput{1}=tElapsed;
%     'svmm': OtherOutput{1}=tElapsed;
%     'hdlm': OtherOutput{1}=tElapsed;
%     'nmf', 'rnmf': OtherOutput{1}=tElapsed;
%     'svdd': OtherOutput{1}=tElapsed;
%     'hkm': OtherOutput{1}=tElapsed;
%     'nc': OtherOutput{1}=tElapsed;
%     'lrc': OtherOutput{1}=tElapsed;
%     'logistic': OtherOutput{1}=tElapsed;
%     'bayesian': OtherOutput{1}=tElapsed;
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

if nargin<4
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

% class statistics
unikCls=unique(trainClass);
numCl=numel(unikCls);
sizeCls=zeros(numCl,1);
for i=1:numCl
   sizeCls(i)=sum(trainClass==unikCls(i)); 
end

% downsample face images
if isfield(option,'ifDownSample')&& option.ifDownSample
    trainSet=downsample(trainSet,option.p);
end
% if tensor data
option.ifTensor=false;
if size(trainSet,3)>1
    option.ifTensor=true;
    option.numR=size(trainSet,1);
    option.numC=size(trainSet,2);
    trainSet=matrizicing(trainSet,3);
    trainSet=trainSet';
end

% number of training samples
[numFeat,numTr]=size(trainSet);

% handle missing value, imputation
if option.ifMissValueImpute
    tfTrain=isnan(trainSet);
    ifMissValTrainSet=any(any(tfTrain));
    if ifMissValTrainSet
        switch option.imputeMethod
            case 'knn'
                trainSet=knnimpute(trainSet);
                model.trainSet=trainSet;
            case 'zero'
                trainSet(tfTrain)=0;
        end
    end
end

% normalization
if logical(option.normalization)
    switch option.normMethod
        case 'mean0std1'
            [trainSet,trainSetMean,trainSetSTD]=normmean0std1(trainSet');
            trainSet=trainSet';
            option.trainSetMean=trainSetMean;
            option.trainSetSTD=trainSetSTD;
        case 'unitl2norm'
            trainSet=normc(trainSet);
    end
end

% classification
switch method
    case 'naiveBayes'
        tic;
        modelBayes = NaiveBayes.fit(trainSet',trainClass);
        model.modelBayes=modelBayes;
        model.method='naiveBayes';
        model.trainOpt=option;
        tElapsed=toc;
        OtherOutput{1}=tElapsed;
    case 'knn'
        optionDefault.k=1;
        option=mergeOption(option,optionDefault);
        tic;
        model.method='knn';
        model.trainOpt=option;
        model.trainSet=trainSet;
        model.trainClass=trainClass;
        tElapsed=toc;
        OtherOutput{1}=tElapsed;
    case 'nc'
        % model selection
        if option.search
            switch option.kernel
                case 'rbf'
                    optGS.classifier='nc';
                    defaultRBFParam=log2(sqrt(numFeat)); % default parameter of rbf
                    optGS.alphaRange=defaultRBFParam-2:0.5:defaultRBFParam+2;
                    optGS.betaRange=0:0;
                    optGS.gammaRange=0:0;
                    optGS.optCl=option;
                    optGS.optCl.search=false;
                    optGS.optCl.normalization=0;
                    if numel(optGS.alphaRange)>1
                    [alphaOptimal,betaOptimal,gammaOptimal,bestAcc,AllAccs]=threeDSearchUniverse(trainSet,trainClass,optGS);
                    option.param=2^alphaOptimal;
                    else
                        option.param=2^optGS.alphaRange;
                    end
                case 'linear'
            end
        end
        tic;
        model=nearestCentroidTrain(trainSet,trainClass,option);
        model.trainOpt=option;
        model.method='nc';
        tElapsed=toc;
        OtherOutput{1}=tElapsed;
%     case 'svm' % libsvm matlab toolbox is required
%         optionDefault.kernel='linear';
%         optionDefault.C=2^0;
%         optionDefault.param=[];
%         option=mergeOption(option,optionDefault);
%         switch option.kernel
%             case 'linear'
%                 % model selection
%                 if option.search % grid search C, sigma
%                     optionGD.classifier='svm';
%                     optionGD.alphaRange=-2:2; % exponent of 2
%                     optionGD.betaRange=0:0; % exponent of 2
%                     optionGD.kfold=3;
%                     optionGD.rerun=1;
%                     optionGD.optCl=option;
%                     optionGD.optCl.normalization=0;
%                     [option.C,option.param]=gridSearchUniverse(trainSet,trainClass,optionGD);
%                 end
%                 if isempty(option.C)
%                    option.trainSetting='-s 0 -t 0 -c 1 -b 1'; 
%                 else
%                    option.trainSetting=['-s 0 -t 0 -c ',num2str(option.C),'  -b 1'];
%                 end
%                 if option.search % line search C
%                     
%                     
%                 end
%             case 'rbf'
%                 % model selection
%                 if option.search % grid search C, sigma
%                     optionGD.classifier='svm';
%                     optionGD.alphaRange=-2:5; % exponent of 2
%                     optionGD.betaRange=3:8; % exponent of 2
%                     optionGD.kfold=3;
%                     optionGD.rerun=1;
%                     optionGD.optCl=option;
%                     optionGD.optCl.normalization=0;
%                     [option.C,option.param]=gridSearchUniverse(trainSet,trainClass,optionGD);
%                 end
%                 if isempty(option.C)
%                    option.trainSetting='-s 0 -t 2 -c 1 -b 1'; 
%                 else
%                     if ~isempty(option.param)
%                         option.trainSetting=['-s 0 -t 2 -c ',num2str(option.C),'  -b 1'];
%                     else
%                         option.trainSetting=['-s 0 -t 2 -g ',num2str(option.param),' -c ',num2str(option.C),'  -b 1'];
%                     end
%                 end
%            otherwise 
%                 error('Please select the correct kernel!');
%         end
%         tic;
%         modelSVM=svmtrain(trainClass,trainSet',option.trainSetting);
%         model.modelSVM=modelSVM;
%         model.trainOpt=option;
%         model.method='svm';
%         tElapsed=toc;
%         OtherOutput{1}=tElapsed;
%     case 'lssvm'
%         tic;
%         testClassPredicted=lssvmClassifier(trainSet,testSet,trainClass,option);
%         tElapsed=toc;
%         OtherOutput{1}=tElapsed;
%     case 'svmLi'     
%         % model selection
%         if option.search
%             switch option.kernel
%                 case 'rbf'
%                     optGS.classifier='svmLi';
%                     optGS.alphaRange=0.05:0.05:1;
%                     optGS.betaRange=-2:1;
%                     optGS.optCl=option;
%                     optGS.optCl.search=false;
%                     optGS.optCl.normalization=0;
%                     [alphaOptimal,betaOptimal,accMax]=gridSearchUniverse(trainSet,trainClass,optGS);
%                     option.C=alphaOptimal;
%                     option.param=2^betaOptimal;
%                 case 'linear'
%                     optGS.classifier='svmLi';
%                     optGS.alphaRange=0.05:0.05:0.95;
%                     optGS.betaRange=0:0;
%                     optGS.optCl=option;
%                     optGS.optCl.search=false;
%                     optGS.optCl.normalization=0;
%                     [alphaOptimal,betaOptimal,accMax]=gridSearchUniverse(trainSet,trainClass,optGS);
%                     option.C=alphaOptimal;
%             end
%         end
%         tic;
%         model=softSVMTrain2(trainSet,trainClass,option);
%        model.trainOpt=option;
%        model.method='svmLi';
%         tElapsed=toc;
% %         perSVM=perform(testClassPredicted,testClass,model.numCl);
%         OtherOutput{1}=tElapsed;
    case 'svmm'
        % model selection
        if option.search
            nuStep=max([0.025,3/(2*mean(sizeCls))]); % step of nu
            maxnu=min(2*min(sizeCls)/numTr-0.05,0.95);
            if nuStep>maxnu
                nuStep=maxnu;
            end
            switch option.kernel
                case 'rbf'
                    optGS.classifier='svmm';
                    defaultRBFParam=log2(sqrt(numFeat)); % default parameter of rbf
                    optGS.alphaRange=nuStep:nuStep:maxnu;%0.025:0.025:maxnu;
                    optGS.betaRange=defaultRBFParam-2:0.5:defaultRBFParam+2;
                    optGS.optCl=option;
                    optGS.optCl.search=false;
                    optGS.optCl.normalization=0;
                    if numel(optGS.alphaRange)>1 || numel(optGS.betaRange)>1
                        [alphaOptimal,betaOptimal,accMax]=gridSearchUniverse(trainSet,trainClass,optGS);
                        option.ctrl=alphaOptimal;
                        option.param=2^betaOptimal;
                    else
                        option.ctrl=optGS.alphaRange;
                        option.param=2^optGS.betaRange;
                    end
                case 'linear'
                    optGS.classifier='svmm';
                    maxnu=min(2*min(sizeCls)/numTr-0.05,0.95);
                    optGS.alphaRange=nuStep:nuStep:maxnu;%0.025:0.025:maxnu;
                    optGS.betaRange=0:0;
                    optGS.optCl=option;
                    optGS.optCl.search=false;
                    optGS.optCl.normalization=0;
                    if numel(optGS.alphaRange)>1
                        [alphaOptimal,betaOptimal,accMax]=gridSearchUniverse(trainSet,trainClass,optGS);
                        option.ctrl=alphaOptimal;
                    else
                        option.ctrl=optGS.alphaRange;
                    end
            end
        end
        tic;
        model=SVMMultiTrain(trainSet,trainClass,option);
       model.trainOpt=option;
       model.method='svmm';
        tElapsed=toc;
%         perSVM=perform(testClassPredicted,testClass,model.numCl);
        OtherOutput{1}=tElapsed;
    case 'svdd'
        % model selection
        if option.search
            nuStep=max([0.025,3/(2*mean(sizeCls))]); % step of nu
            maxnu=0.95;
            switch option.kernel
                case 'rbf'
                    optGS.classifier='svdd';
                    defaultRBFParam=log2(sqrt(numFeat)); % default parameter of rbf
                    optGS.alphaRange=nuStep:nuStep:maxnu;%0.025:0.025:0.95;
                    optGS.betaRange=defaultRBFParam-2:0.5:defaultRBFParam+2;
                    optGS.optCl=option;
                    optGS.optCl.search=false;
                    optGS.optCl.normalization=0;
                    if numel(optGS.alphaRange)>1 || numel(optGS.betaRange)>1
                        [alphaOptimal,betaOptimal,accMax]=gridSearchUniverse(trainSet,trainClass,optGS);
                        option.nu=alphaOptimal;
                        option.param=2^betaOptimal;
                    else
                        option.nu=optGS.alphaRange;
                        option.param=2^optGS.betaRange;
                    end
                case 'linear'
                    optGS.classifier='svdd';
                    optGS.alphaRange=nuStep:nuStep:maxnu;%0.025:0.025:0.95;
                    optGS.betaRange=0:0;
                    optGS.optCl=option;
                    optGS.optCl.search=false;
                    optGS.optCl.normalization=0;
                    if numel(optGS.alphaRange)>1 
                    [alphaOptimal,betaOptimal,accMax]=gridSearchUniverse(trainSet,trainClass,optGS);
                    option.nu=alphaOptimal;
                    else
                        option.nu=optGS.alphaRange;
                    end
            end
        end
        tic;
        model=SVDDTrain(trainSet,trainClass,option);
        model.trainOpt=option;
        model.method='svdd';
        tElapsed=toc;
%         perSVM=perform(testClassPredicted,testClass,model.numCl);
        OtherOutput{1}=tElapsed;
    case 'svddm'
        % model selection
        if option.search
            switch option.kernel
                case 'rbf'
                    optGS.classifier='svddm';
                    optGS.alphaRange=0.05:0.05:1;
                    optGS.betaRange=0.55:0.05:1;
                    optGS.gammaRange=defaultRBFParam-2:0.5:defaultRBFParam+2;
                    optGS.optCl=option;
                    optGS.optCl.search=false;
                    optGS.optCl.normalization=0;
                    [alphaOptimal,betaOptimal,gammaOptimal,bestAcc,AllAccs]=threeDSearchUniverse(trainSet,trainClass,optGS);
                    option.nu=alphaOptimal;
                    option.nus=betaOptimal;
                    option.param=2^gammaOptimal;
                case 'linear'
                    optGS.classifier='svddm';
                    optGS.alphaRange=0.05:0.05:1;
                    optGS.betaRange=0.55:0.05:1;
                    optGS.gammaRange=1:1:1;
                    optGS.optCl=option;
                    optGS.optCl.search=false;
                    optGS.optCl.normalization=0;
                    [alphaOptimal,betaOptimal,gammaOptimal,bestAcc,AllAccs]=threeDSearchUniverse(trainSet,trainClass,optGS);
                    option.nu=alphaOptimal;
                    option.nus=betaOptimal;
            end
        end
        tic;
        model=SVDDMTrain(trainSet,trainClass,option);
        model.trainOpt=option;
        model.method='svddm';
        tElapsed=toc;
%         perSVM=perform(testClassPredicted,testClass,model.numCl);
        OtherOutput{1}=tElapsed;
    case 'svdda'
        % model selection
        if option.search
            switch option.kernel
                case 'rbf'
                    optGS.classifier='svdda';
                    optGS.alphaRange=0.05:0.025:0.95;
                    optGS.betaRange=3:0.5:8;
                    optGS.optCl=option;
                    optGS.optCl.search=false;
                    optGS.optCl.normalization=0;
                    [alphaOptimal,betaOptimal,accMax]=gridSearchUniverse(trainSet,trainClass,optGS);
                    option.fracsv=alphaOptimal;
                    option.param=2^betaOptimal;
                case 'linear'
                    optGS.classifier='svdda';
                    optGS.alphaRange=0.05:0.05:0.95;
                    optGS.betaRange=0:0;
                    optGS.optCl=option;
                    optGS.optCl.search=false;
                    optGS.optCl.normalization=0;
                    [alphaOptimal,betaOptimal,accMax]=gridSearchUniverse(trainSet,trainClass,optGS);
                    option.fracsv=alphaOptimal;
            end
        end
        tic;
        model=SVDDAdaptiveTrain(trainSet,trainClass,option);
        model.trainOpt=option;
        model.method='svdda';
        tElapsed=toc;
%         perSVM=perform(testClassPredicted,testClass,model.numCl);
        OtherOutput{1}=tElapsed;
    case 'lsvm' % local svm
        tic;
        model.method='lsvm';
        model.trainOpt=option;
        model.trainSet=trainSet;
        model.trainClass=trainClass;
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
        optionDefault=[];
        option=mergeOption(option,optionDefault);
        % model selection % search parameters for svmLi, nnls, hdlm
        if option.search
            switch option.classMethod
                case {'svmLi','svm'}
                    if isfield(option.optCl,'kernel') && strcmp(option.optCl.kernel,'linear')
                    optGS.alphaRange=-3:5;
                    optGS.betaRange=0:0;
                    end
                    if isfield(option.optCl,'kernel') && strcmp(option.optCl.kernel,'rbf')
                    optGS.alphaRange=2:5;
                    optGS.betaRange=6:10;
                    end
                case {'nnls','hdlm'}
                    if isfield(option.optCl,'kernel') && strcmp(option.optCl.kernel,'rbf')
                        optGS.alphaRange=0:0;
                        optGS.betaRange=4:9;
                    end
            end
            optGS.classifier='hkm';
            optGS.optCl=option;
            optGS.optCl.search=false;
            optGS.optCl.normalization=0;
            tic;
            [alphaOptimal,betaOptimal,accMax]=gridSearchUniverse(trainSet,trainClass,optGS);
            option.optCl.C=2^alphaOptimal;
            option.optCl.param=2^betaOptimal;
        end
        model=hierarchialModelTrain(trainSet,trainClass,[],option);
        model.trainOpt=option;
        model.method='hkm';
        tElapsed=toc;
        OtherOutput{1}=tElapsed;
%     case 'bayesian'
%         optionDefault.type='diaglinear'; % 'linear','diaglinear','quadratic','diagquadratic','mahalanobis'
%         option=mergeOption(option,optionDefault);
%         tic;
%         [testClassPredicted,trainError,testPosterior]=classify(testSet',trainSet',trainClass,option.type);
%         tElapsed=toc;
%         OtherOutput{1}=tElapsed;
%         OtherOutput{2}=trainError;
%         OtherOutput{3}=testPosterior;
    case 'nnls'
        if isfield(option,'kernel') && strcmp(option.kernel,'ds')
            if isfield(option,'param') && isempty(option.param) || ~isfield(option,'param')
                option.param(1)=1;
                option.param(2)=5;
            end
            option.param=[option.numR;option.numC;option.param(1);option.param(2)];
        end
        % model selection
        if option.search
            optGS.classifier='nnls';
            optGS.alphaRange=0:0;
            optGS.betaRange=3:8;
            [alphaOptimal,betaOptimal,accMax]=gridSearchUniverse(trainSet,trainClass,optGS);
        option.param=2^betaOptimal;
        end
        tic;
        model.trainOpt=option;
        model.method='nnls';
        model.trainSet=trainSet;
        model.trainClass=trainClass;
        tElapsed=toc;
        OtherOutput{1}=tElapsed;
        
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
        optionDefault.numRandom=100;
        option=mergeOption(option,optionDefault);
        tic;
        if strcmp(option.kernel,'ds')
            rank=NaN;
            lambda=NaN;
            option.param=[option.numR;option.numC;rank;lambda];
        end
        model.trainOpt=option;
        model.method='bootstrapnnls';
        model.trainSet=trainSet;
        model.trainClass=trainClass;
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
        optionDefault.facts=numel(unique(trainClass));
        option=mergeOption(option,optionDefault);
        tic;
        model.trainOpt=option;
        model.method='nmf';
        model.trainSet=trainSet;
        model.trainClass=trainClass;
        tElapsed=toc;
        OtherOutput{1}=tElapsed;
    case 'rnmf'
        optionDefault.facts=7;%numel(unique(trainClass));
        option=mergeOption(option,optionDefault);
        tic;
        model.trainOpt=option;
        model.method='rnmf';
        model.trainSet=trainSet;
        model.trainClass=trainClass;
        tElapsed=toc;
        OtherOutput{1}=tElapsed;
%     case 'src'
%         tic;
%         [testClassPredicted,sparsity]=src(trainSet,trainClass,testSet,option);
%         tElapsed=toc;
%         OtherOutput{1}=tElapsed;
%         OtherOutput{2}=sparsity;
%     case 'src2'
%         if isfield(option,'kernel') && strcmp(option.kernel,'ds')
%             if isfield(option,'param') && isempty(option.param) || ~isfield(option,'param')
%                 option.param(1)=1;
%                 option.param(2)=5;
%             end
%             option.param=[option.numR;option.numC;option.param(1);option.param(2)];
%         end
%         if option.search
%             optGS.classifier='src2';
%             optGS.alphaRange=-12:-6;
%             optGS.betaRange=0:0;
%             [alphaOptimal,betaOptimal,accMax]=gridSearchUniverse(trainSet,trainClass,optGS);
%             option.lambda=2^alphaOptimal;
%             option.param=2^betaOptimal;
%         end
%         tic;
%         model.trainOpt=option;
%         model.method='src2';
%         model.trainSet=trainSet;
%         model.trainClass=trainClass;
%         tElapsed=toc;
%         OtherOutput{1}=tElapsed;
    case 'ksrsc'
        if isfield(option,'kernel') && strcmp(option.kernel,'ds')
            if isfield(option,'param') && isempty(option.param) || ~isfield(option,'param')
                option.param(1)=1;
                option.param(2)=5;
            end
            option.param=[option.numR;option.numC;option.param(1);option.param(2)];
        end
        if option.search
            optGS.classifier='ksrsc';
            optGS.alphaRange=-12:-6;
            optGS.betaRange=0:0;
            [alphaOptimal,betaOptimal,accMax]=gridSearchUniverse(trainSet,trainClass,optGS);
            option.lambda=2^alphaOptimal;
            option.param=2^betaOptimal;
        end
        tic;
        model.trainOpt=option;
        model.method='ksrsc';
        model.trainSet=trainSet;
        model.trainClass=trainClass;
        tElapsed=toc;
        OtherOutput{1}=tElapsed;
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
%     case 'msrc'
%         tic;
%         optionDefault.ifModelSelection=false;
%         optionDefault.metaSampleMethod='svd';
%         option=mergeOption(option,optionDefault);
%         [trainSet,trainClass]=computeMetaSample(trainSet,trainClass,option);
%         if option.ifModelSelection
%             option.normalization=0;
%             [param,maxAcc]=lineSearchSRC2(trainSet,trainClass,option);
%             option.lambda=param;
%         end
%         model.method='msrc';
%         model.metaSample=trainSet;
%         model.metaClass=trainClass;
%         model.trainOpt=option;
%         tElapsed=toc;
%         OtherOutput{1}=tElapsed;
%         OtherOutput{2}=0;
%     case 'mnnls'
%         tic;
%         optionDefault.ifModelSelection=false;
%         optionDefault.metaSampleMethod='nmf';
%         option=mergeOption(option,optionDefault);
%         [trainSet,trainClass]=computeMetaSample(trainSet,trainClass,option);
%         [testClassPredicted,sparsity]=nnlsClassifier(trainSet,trainClass,testSet,testClass,option);
% %         methodCl='svm'; % SVM does not work
% %         optionsvm.normalization=0;
% %         optionsvm.trainSetting='-t 2 -c 1 -b 1';
% %         optionsvm.testSetting='-b 1';
% %         optionsvm.search=false;
% %         [testClassPredicted,classPerform,OtherOutput]=classification(trainSet,trainClass,testSet,testClass,methodCl,optionsvm);
% %         sparsity=[];
%         tElapsed=toc;
%         OtherOutput{1}=tElapsed;
%         OtherOutput{2}=sparsity;
    case 'subdic'
        tic;
        optionDefault.ifModelSelection=false;
        optionDefault.metaSampleMethod='vsmf';
        option=mergeOption(option,optionDefault);
        [trainSet,trainClass,optionReturned]=computeMetaSample(trainSet,trainClass,option);
%         if option.ifModelSelection
%             option.normalization=0;
%             [param,maxAcc]=lineSearchSRC2(trainSet,trainClass,option);
%             option.lambda=param;
%         end
        model.method='ksrsc';
        model.trainSet=trainSet;
        model.trainClass=trainClass;
        model.trainOpt=optionReturned;
        tElapsed=toc;
        OtherOutput{1}=tElapsed;
        OtherOutput{2}=0;
    case 'lrc'
        % model selection
        if option.search
            switch option.kernel
                case 'rbf'
                    optGS.classifier='lrc';
                    defaultRBFParam=log2(sqrt(numFeat)); % default parameter of rbf
                    optGS.alphaRange=defaultRBFParam-2:0.5:defaultRBFParam+2;
                    optGS.betaRange=0:0;
                    optGS.gammaRange=0:0;
                    optGS.optCl=option;
                    optGS.optCl.search=false;
                    optGS.optCl.normalization=0;
                    if numel(optGS.alphaRange)>1
                    [alphaOptimal,betaOptimal,gammaOptimal,bestAcc,AllAccs]=threeDSearchUniverse(trainSet,trainClass,optGS);
                    option.param=2^alphaOptimal;
                    else
                        option.param=2^optGS.alphaRange;
                    end
                case 'linear'
            end
        end
        tic;
        model.trainOpt=option;
        model.method='lrc';
        model.trainSet=trainSet;
        model.trainClass=trainClass;
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
        % model selection
        if option.search && strcmp(option.kernel,'rbf')
            optGS.classifier='hdlm';
            optGS.alphaRange=0:0;
            optGS.betaRange=4:9;
            optGS.optCl=option;
            optGS.optCl.search=false;
            optGS.optCl.normalization=0;
            [alphaOptimal,betaOptimal,accMax]=gridSearchUniverse(trainSet,trainClass,optGS);
            option.param=2^betaOptimal;
        end
        if option.search && strcmp(option.kernel,'sigmoid')
            optGS.classifier='hdlm';
            optGS.alphaRange=-3:3;
            optGS.betaRange=-3:3;
            optGS.optCl=option;
            optGS.optCl.search=false;
            optGS.optCl.normalization=0;
            [alphaOptimal,betaOptimal,accMax]=gridSearchUniverse(trainSet,trainClass,optGS);
            option.param=[alphaOptimal,betaOptimal];
        end
        model=khdlmTrain(trainSet,trainClass,option);
        model.trainOpt=option;
        model.method='hdlm';
        tElapsed=toc;
        OtherOutput{1}=tElapsed;
    case 'elm'
        tic;
        Elm_Type=1;
        optionDefault.TimesofHiddenNeurons=3;
        optionDefault.ActivationFunction='sig';
        option=mergeOption(option,optionDefault);
        option.NumberofHiddenNeurons=option.TimesofHiddenNeurons * size(trainSet,2);
        %         if option.NumberofHiddenNeurons==0
        %             option.NumberofHiddenNeurons=size(trainSet,2);
        %         end
        model.method='elm';
        model.trainSet=trainSet;
        model.trainClass=trainClass;
        model.trainOpt=option;
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
        optionDefault=[];
        option=mergeOption(option,optionDefault);
        model=logisticRegressTrain(trainSet,trainClass,option);
        model.trainOpt=option;
        model.method='logistic';
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