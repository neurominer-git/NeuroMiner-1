function [testClassPredicted, sparsity,OtherOuptput] = SRC2(trainSet, trainClass, testSet, option)
% SRC2 for overcomplete data
% trainSet: matrix, each column is a training sample
% trainClass: column vector, the class labels for training samples
% testSet: matrix, the test samples
% option, struct, with fields:
% option.lambda: scalar, the parameter to optimization algorithm l1_ls, the
% default is 0.1.
% option.predicter: string, the rule to interpret the sparse code, it can be
% 'subspace' (default),'max','kvote'.
% option.method: string, the optimization method. It can be 'activeSet',
% 'interiorPoint', 'proximal', and 'smo'.
% testClassPredicted: column vectors, the predicted class labels of
% testing samples.
% sparsity: scalar, the sparsity of the sparse coefficient matrix.
% Yifeng Li
% note: each sample has to be normalized to unit l2 norm

% % normalization to length 1
% trainSet=normc(trainSet);
% testSet=normc(testSet);

% optionDefault.p=4;
optionDefault.lambda=0.1;
optionDefault.kernel='linear';
optionDefault.param=0;
optionDefault.predicter='subspace';
optionDefault.knn=numel(trainClass);
optionDefault.method='activeSet';
option=mergeOption(option,optionDefault);
% trainSet=downsample(trainSet,option.p);
% testSet=downsample(testSet,option.p);

% training step, obtain sparse coefficients in columns of Y
Y=zeros(size(trainSet,2),size(testSet,2));
switch option.method
    case 'activeSet'
        % active set algorithm
        H=computeKernelMatrix(trainSet,trainSet,option);
        C=computeKernelMatrix(trainSet,testSet,option);
        [Y,~]=l1QPActiveSet(H,-C,option.lambda);
    case 'interiorPoint'
        % interior point-method
%         AtA=computeKernelMatrix(trainSet,trainSet,option);
%         AtB=computeKernelMatrix(trainSet,testSet,option);
%         BtB=computeKernelMatrix(testSet,testSet,option);
%         BtB=diag(BtB);
%         Y=l1LSKernelBatchDL(AtA,AtB,BtB,option);
        for i=1:size(testSet,2)
            Y(:,i)= l1_ls(trainSet, testSet(:,i), option.lambda); %http://www.stanford.edu/~body/l1_ls
        end
    case 'proximal'
        % proximal method
        AtA=computeKernelMatrix(trainSet,trainSet,option);
        AtB=computeKernelMatrix(trainSet,testSet,option);
        BtB=computeKernelMatrix(testSet,testSet,option);
        BtB=diag(BtB);
        Y=l1LSProximal(AtA,AtB,BtB,option);
    case 'smo'
        H=computeKernelMatrix(trainSet,trainSet,option);
        C=computeKernelMatrix(trainSet,testSet,option);
        Y=l1QPSMOMulti(H,-C,option.lambda);
    otherwise 
        error('choose correct method for l1ls');
end
% calculate sparsity
sparsity=sum(sum(abs(Y)<=0.0001))/(size(Y,1)*size(Y,2));

% predict step
switch option.predicter
    case  'max'
        [val,ind]=max(Y,[],1);
        testClassPredicted=trainClass(ind);
        OtherOuptput=[];
    case 'knn'
        testClassPredicted=knnrule(Y,trainClass,option.knn);
        OtherOuptput=[];
    case {'subspace','ns'}
        [testClassPredicted,residuals]=subspace(Y,testSet,trainSet,trainClass);
        OtherOuptput=residuals;
end
end
