function [trainExtr,out,time]=sr_featureExtractionTrain(trainSet,trainClass,feMethod,opt)
% Extract features from training set
% Usuage:
% [trainExtr,out]=featureExtractionTrain(trainSet,trainClass,feMethod)
% [trainExtr,out]=featureExtractionTrain(trainSet,trainClass,feMethod,opt)
% trainSet, matrix, the training set with samples in columns and features in rows.
% trainClass: column vector of numbers or string, the class labels of the traning set.
% feMethod: string, feature extraction methods. It could be
%     'none': there is no feature extraction, return the input trainSet;
%     'ksrdl': kernel sparse representation/dictionary learning based method;
%     'vsmf': versatile sparse matrix factorization;
%     'nmf': the standard NMF;
%     'sparsenmf': sparse-NMF;
%     'orthnmf': orthogonal-NMF;
%     'seminmf': semi-NMF;
%     'convexnmf': convexc-NMF;
%     'knmf-dc': knmf-dc;
%     'knmf-cv': knmf-cv, kernel convex NMF;
%     'knmf-nnls': knmf-nnls, kernel semi-NMF using NNLS algorithm;
%     'knmf-ur': knmf-ur, kernel semi-NMF using update-rule algorithm;
%     'pca': PCA.
% opt: struct, the options to configure the feature extraction.
% opt.facts: scalar, the number of new features to extract.
% opt.option: options for specific algorithm. If feMethod is
%     'non': opt.option is invalid;
%     'ksrdl': kernel sparse representation/dictionary learning based method;
%     'vsmf': versatile sparse matrix factorization;
%     'nmf': type "help nmfnnls" for more information;
%     'sparsenmf': type "help sparsenmfnnls" for more information;
%     'orthnmf': type "help orthnmfrule" for more information;
%     'seminmf': type "help seminmfnnls" for more information;
%     'convexnmf': type "help convexnmfrule" for more information;
%     'knmf-dc': type "help kernelnmfdecom" for more information;
%     'knmf-cv': type "help kernelconvexnmf" for more information;
%     'knmf-nnls': type "help kernelseminmfnnls" for more information;
%     'knmf-ur': type "help kernelseminmfrule" for more information;
%     'pca': opt.option is invalid.
% trainExtr: matrix, the training data in the feature space.
% out: struct, include the feature extraction model information.
%     out.feMethod: the same as the input argument "feMethod";
%     out.facts: the same as opt.fact;
%     out.factors: column vector of cell, the factor matrices produced by specific algorithm;
%     out.optionTr: struct, = opt.option.
% Reference:
% [1]\bibitem{bibm10}
%     Y. Li and A. Ngom,
%     ``Non-negative matrix and tensor factorization based classification of clinical microarray gene expression data,''
%     {\it IEEE International Conference on Bioinformatics \& Biomedicine},
%     2010, pp.438-443.
% [2]\bibitem{cibcb2012}
%     Y. Li and A. Ngom,
%     ``A New Kernel Non-Negative Matrix Factorization and Its Application in Microarray Data Analysis,''
%     {\it CIBCB},
%     submited, 2012.
% Contact Information:
% Yifeng Li
% University of Windsor
% li11112c@uwindsor.ca; yifeng.li.cn@gmail.com
% May 22, 2011
if nargin<4
    opt=[];
end
optDefault.facts=numel(unique(trainClass));
optDefault.normalization=0;
optDefault.normMethod='';
optDefault.option=[];
opt=mergeOption(opt,optDefault);

out.normalization=opt.normalization;
out.normMethod=opt.normMethod;
% normalization
if logical(opt.normalization)
    switch opt.normMethod
        case 'mean0std1'
            [trainSet,trainSetMean,trainSetSTD]=normmean0std1(trainSet');
            trainSet=trainSet';
            out.trainSetMean=trainSetMean;
            out.trainSetSTD=trainSetSTD;
        case 'unitl2norm'
            trainSet=normc(trainSet);
    end
end

tStart=tic;
switch feMethod
    case 'none'
        trainExtr=trainSet;
        out.feMethod=feMethod;
    case 'nmf'
        %         if ~isfield(opt,'facts')
        %            opt.facts=numel(unique(trainClass));
        %         end
        %         optionDefault.option=[];
        %         opt=mergeOption(opt,optionDefault);
        [A,trainExtr]=nmfnnls(trainSet,opt.facts,opt.option);
        out.feMethod=feMethod;
        out.factors={A;trainExtr};
        out.facts=opt.facts;
        out.optionTr=opt.option;
    case 'sparsenmf'
        %         if ~isfield(opt,'facts')
        %            opt.facts=numel(unique(trainClass));
        %         end
        %         optionDefault.option=[];
        %         opt=mergeOption(opt,optionDefault);
        [A,trainExtr]=sparsenmfnnls(trainSet,opt.facts,opt.option);
        out.feMethod=feMethod;
        out.factors={A;trainExtr};
        out.facts=opt.facts;
        out.optionTr=opt.option;
    case 'orthnmf'
        %         if ~isfield(opt,'facts')
        %            opt.facts=numel(unique(trainClass));
        %         end
        %         optionDefault.option=[];
        %         opt=mergeOption(opt,optionDefault);
        [A,S,trainExtr]=orthnmfrule(trainSet,opt.facts,opt.option);
        out.feMethod=feMethod;
        out.factors={A;S;trainExtr};
        out.facts=opt.facts;
        out.optionTr=opt.option;
    case 'knmf-nnls'
        %         if ~isfield(opt,'facts')
        %             opt.facts=numel(unique(trainClass));
        %         end
        %         optionDefault.option=[];
        %         opt=mergeOption(opt,optionDefault);
        [AtA,trainExtr,~,~,~]=kernelseminmfnnls(trainSet,opt.facts,opt.option);
        out.feMethod=feMethod;
        out.factors={AtA;trainExtr;trainSet};
        out.facts=opt.facts;
        out.optionTr=opt.option;
    case 'knmf-l1nnls'
        %         if ~isfield(opt,'facts')
        %             opt.facts=numel(unique(trainClass));
        %         end
        %         optionDefault.option=[];
        %         opt=mergeOption(opt,optionDefault);
        [AtA,trainExtr,~,~,~]=kernelsparseseminmfnnls(trainSet,opt.facts,opt.option);
        out.feMethod=feMethod;
        out.factors={AtA;trainExtr;trainSet};
        out.facts=opt.facts;
        out.optionTr=opt.option;
    case 'knmf-l1nnqp'
        %         if ~isfield(opt,'facts')
        %             opt.facts=numel(unique(trainClass));
        %         end
        %         optionDefault.option=[];
        %         opt=mergeOption(opt,optionDefault);
        [AtA,trainExtr,~,~,~]=kernelSparseNMFNNQP(trainSet,opt.facts,opt.option);
        out.feMethod=feMethod;
        out.factors={AtA;trainExtr;trainSet};
        out.facts=opt.facts;
        out.optionTr=opt.option;
    case 'knmf-ur'
        [AtA,trainExtr,~,~,~]=kernelseminmfrule(trainSet,opt.facts,opt.option);
        out.feMethod=feMethod;
        out.factors={AtA;trainExtr;trainSet};
        out.facts=opt.facts;
        out.optionTr=opt.option;
    case 'knmf-dc'
        [A,trainExtr,~,~,~]=kernelnmfdecom(trainSet,opt.facts,opt.option);
        out.feMethod=feMethod;
        out.factors={A;trainExtr};
        out.facts=opt.facts;
        out.optionTr=opt.option;
    case 'knmf-cv'
        [XtX,W,trainExtr,~,~,~]=kernelconvexnmf(trainSet,opt.facts,opt.option);
        AtA=W'*XtX*W;
        out.feMethod=feMethod;
        out.factors={AtA;W;trainExtr;trainSet};
        out.facts=opt.facts;
        out.optionTr=opt.option;
    case {'ksr','ksrdl'}
        % dictionary learning
        trainTrainU=computeKernelMatrix(trainSet,trainSet,opt.option);
        if strcmp(opt.option.kernel,'ds')
           trainTrain=normalizeKernelMatrix(trainTrainU);
           [AtA,trainExtr,~,pinvY,~,~,~]=KSRDL([],trainTrain,opt.facts,opt.option);
        else
            [AtA,trainExtr,~,pinvY,~,~,~]=KSRDL([],trainTrainU,opt.facts,opt.option);
        end
        out.feMethod=feMethod;
        out.factors={AtA;trainExtr;trainSet;trainTrainU};
        out.pinvY=pinvY;
        out.facts=opt.facts;
        out.optionTr=opt.option;
%     case 'ksr-l1ls'
%         opt.option.SRMethod='l1ls';
%         % dictionary learning
%         trainTrain=computeKernelMatrix(trainSet,trainSet,opt.option);
%         [AtA,trainExtr,~,~,~]=KSRDL(trainTrain,opt.facts,opt.option);
%         out.feMethod=feMethod;
%         out.factors={AtA;trainExtr;trainSet};
%         out.facts=opt.facts;
%         out.optionTr=opt.option;
%     case 'ksr-l1qpAs'
%         opt.option.SRMethod='l1qpAs';
%         % dictionary learning
%         trainTrain=computeKernelMatrix(trainSet,trainSet,opt.option);
%         [AtA,trainExtr,~,~,~]=KSRDL(trainTrain,opt.facts,opt.option);
%         out.feMethod=feMethod;
%         out.factors={AtA;trainExtr;trainSet};
%         out.facts=opt.facts;
%         out.optionTr=opt.option;
%     case 'ksr-proxl1ls'
%         opt.option.SRMethod='proxl1ls';
%         % dictionary learning
%         trainTrain=computeKernelMatrix(trainSet,trainSet,opt.option);
%         [AtA,trainExtr,~,~,~]=KSRDL(trainTrain,opt.facts,opt.option);
%         out.feMethod=feMethod;
%         out.factors={AtA;trainExtr;trainSet};
%         out.facts=opt.facts;
%         out.optionTr=opt.option;
%     case 'ksrdr-l1ls'
%         opt.option.SRMethod='l1ls';
%         % dictionary learning
%         trainTrain=computeKernelMatrix(trainSet,trainSet,opt.option);
%         % normalize trainTrain
%         %         trainTrain=normc(trainTrain);
%         optK.kernel='linear';
%         trainTrain2=computeKernelMatrix(trainTrain,trainTrain,optK);
%         [AtA,trainExtr,~,~,~]=KSRDLDR(trainTrain2,opt.facts,opt.option);
%         out.feMethod=feMethod;
%         out.factors={AtA;trainExtr;trainSet};
%         out.facts=opt.facts;
%         out.optionTr=opt.option;
%     case 'ksr-nnls'
%         opt.option.SRMethod='nnls';
%         % dictionary learning
%         trainTrain=computeKernelMatrix(trainSet,trainSet,opt.option);
%         [AtA,trainExtr,~,~,~]=KSRDL(trainTrain,opt.facts,opt.option);
%         out.feMethod=feMethod;
%         out.factors={AtA;trainExtr;trainSet};
%         out.facts=opt.facts;
%         out.optionTr=opt.option;
%     case 'ksr-l1nnls'
%         opt.option.SRMethod='l1nnls';
%         % dictionary learning
%         trainTrain=computeKernelMatrix(trainSet,trainSet,opt.option);
%         [AtA,trainExtr,~,~,~]=KSRDL(trainTrain,opt.facts,opt.option);
%         out.feMethod=feMethod;
%         out.factors={AtA;trainExtr;trainSet};
%         out.facts=opt.facts;
%         out.optionTr=opt.option;
    case 'seminmf'
        %         if ~isfield(opt,'facts')
        %            opt.facts=numel(unique(trainClass));
        %         end
        %         optionDefault.option=[];
        %         opt=mergeOption(opt,optionDefault);
        [A,trainExtr]=seminmfnnls(trainSet,opt.facts,opt.option); % seminmfrule
        out.feMethod=feMethod;
        out.factors={A;trainExtr};
        out.facts=opt.facts;
        out.optionTr=opt.option;
    case 'convexnmf'
        %         if ~isfield(opt,'facts')
        %            opt.facts=numel(unique(trainClass));
        %         end
        %         optionDefault.option=[];
        %         opt=mergeOption(opt,optionDefault);
        [A,trainExtr]=convexnmfrule(trainSet,opt.facts,opt.option);
        out.feMethod=feMethod;
        out.factors={trainSet;A;trainExtr};
        out.facts=opt.facts;
        out.optionTr=opt.option;
    case 'transductive-nmf'
        numTr=size(trainSet,2);
%         numTe=size(testSet,2);
        [AtA,trainTestExtr,~,~,~]=kernelseminmfnnls([trainSet,testSet],opt.facts,opt.option);
        out.feMethod=feMethod;
        trainExtr=trainTestExtr(:,1:numTr);
        testExtr=trainTestExtr(:,numTr+1:end);
        out.factors={AtA;trainExtr;trainSet;testExtr};
        out.facts=opt.facts;
        out.optionTr=opt.option;
    case 'vsmf'
        [A,trainExtr,AtA]=vsmf(trainSet,opt.facts,opt.option);
        out.feMethod=feMethod;
        out.factors={AtA;A;trainExtr;trainSet};
        out.facts=size(AtA,1);
        out.optionTr=opt.option;
    case 'pca'
        %         optionDefault.option=[];
        %         opt=mergeOption(opt,optionDefault);
        A = princomp(trainSet','econ');
        A=A(:,1:opt.facts);
        trainExtr=A'* trainSet;
        out.feMethod=feMethod;
        out.factors={A;trainExtr};
        out.facts=opt.facts;
        out.optionTr=opt.option;
    otherwise
        error('Please find the correct feature extraction method or define your own.');
        
end
time=toc(tStart);
out.trainSet=trainSet;
out.trainClass=trainClass;
end