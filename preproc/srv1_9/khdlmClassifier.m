function [testClassPredicted,beta]=khdlmClassifier(trainSet,testSet,trainClass,testClass,option)
% High Dimensional Linear Machine with multiple output nodes
% option: struct:
% option.kernel: string, kernel function. It could be
%     'linear': linear kernel, linear(x,y)=x'*y;
%     'polynomial': polynomial kernel, polynomial(x,y)=(Gamma.*(x'*y)+ Coefficient).^Degree;
%     'rbf': rdial basis function kernel (default), rbf(x,y)=exp((-(1/sigma^2)).*(|x-y|.^2);
%     'sigmoid': sigmoid kernel, sigmoid(x,y)=tanh(alpha*(x'*y) + beta).
%     yourKernelFunctionName: you can define and use you own kernel function 
% option.param: parameters of specified kernel:
%      if option.kernel='linear': option.param=[];
%      if option.kernel='polynomial': option.param is [Gamma;Coefficient;Degree], the default is [1;0;2];
%      if option.kernel='rbf': option.param is sigma, the default is 1/size(X,1);
%      if option.kernel='sigmoid': option.param is [alpha;beta], the default is [1;0];
%      if option.kernel use your own kernel function, option.param is the parameter of your kernel.

%August 28,2011

optionDefault.kernel='linear';
optionDefault.param=[];
if nargin<5
   option=optionDefault;
else
    option=mergeOption(option,optionDefault);
end
switch option.kernel
    case 'rbf'
        if isempty(option.param)
            option.param=1/size(trainSet,1);
        end
%     sigma=param(1);
%     kfun= @kernelRBF; % my rbf kernel
    kfun=@kernelRBF; % fast rbf kernel from official matla
    case 'polynomial'
        if isempty(option.param)
            option.param=[1;0;2];
        end
%     Gamma=param(1);
%     Coefficient=param(2);
%     Degree=param(3);
    kfun= @kernelPoly;
    case 'linear'
    kfun= @kernelLinear;
    case 'sigmoid'
        if isempty(option.param)
            option.param=[1;0];
        end
%         alpha=param(1);
%         beta=param(2);
        kfun=@kernelSigmoid;
    otherwise
        eval(['kfunc=@',option.kernel,';']);
end

numTr=size(trainSet,2);
numTe=size(testSet,2);
unikClass=unique([trainClass;testClass]);
unikClass=sort(unikClass);
numUnik=numel(unikClass);
trainClassM=-ones(numTr,numUnik);
trainSetM=trainSet;
TrPointer=1; % the next position to insert
for i=1:numUnik
    curTrInd=(trainClass==unikClass(i));
    numCurTr=sum(curTrInd);
    trainSetM(:,TrPointer:TrPointer+numCurTr-1)=trainSet(:,curTrInd);
    trainClassM(TrPointer:TrPointer+numCurTr-1,i)=1;
    TrPointer=TrPointer+numCurTr;
end

trainSetM=[ones(1,numTr);trainSetM];
testSetM=[ones(1,numTe);testSet];
clear('trainSet','testSet');
% compute kernel matrix
% K=feval(kfun,trainSetM,trainSetM,option.param); % kernel matrix
K=computeKernelMatrix(trainSetM,trainSetM,option);
% Kte=feval(kfun,trainSetM,testSetM,option.param); % kernel matrix
Kte=computeKernelMatrix(trainSetM,testSetM,option);
% weights
% testClassPredictedM=trainClassM'*(K\Kte);
testClassPredictedM=trainClassM'*invsvd(K)*Kte;
testClassPredicted=nan(numTe,1);
for i=1:size(testClassPredicted)
     [valp,indp]=max(testClassPredictedM(:,i));
     testClassPredicted(i)=unikClass(indp);
%      [val,ind]=max(testClassM(i,:));
%      if indp==ind
%          testClassPredicted(i)=unikClass(ind);
%      end
end
% normal vector
beta=pinv(trainSetM')*trainClassM(:,1);
% beta=trainSetM*invsvd(K)*trainClassM(:,1);
end

