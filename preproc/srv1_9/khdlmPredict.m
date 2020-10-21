function testClassPredicted=khdlmPredict(model,testSet)
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

numTe=size(testSet,2);
testSetM=[ones(1,numTe);testSet];
clear('testSet');
% compute kernel matrix;
% Kte=feval(kfun,trainSetM,testSetM,option.param); % kernel matrix
Kte=computeKernelMatrix(model.trainSetM,testSetM,model.option);
% weights
% testClassPredictedM=trainClassM'*(K\Kte);
testClassPredictedM=model.trainClassMtKInv * Kte;
testClassPredicted=nan(numTe,1);
for i=1:size(testClassPredicted)
     [valp,indp]=max(testClassPredictedM(:,i));
     testClassPredicted(i)=model.unikClass(indp);
end
end

