function [eigvector, elapse, Y] = SDA_MR(gnd,fea,semiSplit,options,RegularizerOptions)
% SDA: Semi-supervised Discriminant Analysis (Manifold Regularization Way)
%       NOTE: Ref [2][3] shows that LDA can be formulated as a regression
%             problem. Thus, manifold regularization [3][4] can naturally
%             be used to get the semi-supervised discriminant analysis
%             algorithm.  
%
%       [eigvector, eigvalue, elapse] = SDA_MR(gnd,fea,semiSplit,KernelOptions,RegularizerOptions)
% 
%             Input:
%               gnd     - Label vector.  
%               fea     - data matrix. Each row vector of fea is a data point. 
%
%             semiSplit - fea(semiSplit,:) is the labeled data matrix. 
%                       - fea(~semiSplit,:) is the unlabeled data matrix. 
%
%
%               options - Struct value in Matlab. The fields in options
%                         that can be set:
%
%                     WOptions     Please see ConstructW.m for detailed options.  
%                       or
%                      W           You can construct the W outside.
%
%           LaplacianNorm = 0 | 1  (0 for un-normalized and 1 for
%                                   normalized graph laplacian) 
%                                   Default: 0
%             LaplacianDegree      power of the graph Laplacian to use as
%                                  the graph regularizer
%                                   Default: 1
%
%                   ReguBeta         Paramter for manifold regularizer
%                                    Default 0.1 
%
%                         
%           KernelType  -  Choices are:
%               'Gaussian'      - e^{-(|x-y|^2)/2t^2}
%               'Polynomial'    - (x'*y)^d
%               'PolyPlus'      - (x'*y+1)^d
%               'Linear'        -  x'*y
%
%               t       -  parameter for Gaussian
%               d       -  parameter for Poly
%
%             Output:
%               eigvector - Each column is an embedding function, for a new
%                           data point (row vector) x,  
%                                   y = x*eigvector (KernelType = 'Linear')
%                                   y = K(x,:)*eigvector (KernelType = other)
%                                       K(x,:) = [K(x1,x),K(x2,x),...K(xm,x)]
%                           will be the embedding result of x.
%               elapse    - Time spent on different steps 
%               Y         - Embedding results
% 
%
%    Examples:
%
%       
% 
%
% See also LPP, LGE
%
%Reference:
%
%   [1] Deng Cai, Xiaofei He and Jiawei Han, "Semi-Supervised Discriminant
%   Analysis ", IEEE International Conference on Computer Vision (ICCV),
%   Rio de Janeiro, Brazil, Oct. 2007.   
%
%   [2] Deng Cai, Xiaofei He, Jiawei Han, "Efficient Kernel Discriminant
%   Analysis via Spectral Regression", Proc. 2007 Int. Conf. on Data
%   Mining (ICDM'07), Omaha, NE, Oct. 2007. 
%
%   [3] Deng Cai, Xiaofei He and Jiawei Han, "SRDA: An Efficient Algorithm for
%   Large Scale Discriminant Analysis" IEEE Transactions on Knowledge and
%   Data Engineering, January 2008.  
%
%   [4] M. Belkin, P. Niyogi, V. Sindhwani, "Manifold Regularization: A
%   Geometric Framework for Learning from Labeled and Unlabeled Examples",
%   Journal of Machine Learning Research 7(Nov):2399--2434, 2006 
%
%   [5] V. Sindhwani, P. Niyogi, M. Belkin, "Beyond  the  Point  Cloud:  from
%   Transductive  to  Semi-supervised  Learning", ICML 2005.
%
%   version 2.0 --November/2007 
%   version 1.0 --May/2007 
%
%   Written by Deng Cai (dengcai2 AT cs.uiuc.edu)


ReguBeta = 0.1;
if isfield(options,'ReguBeta') 
    ReguBeta = options.ReguBeta;
end


[nSmp,nFea] = size(fea);
nSmpLabel = sum(semiSplit);
nSmpUnlabel = sum(~semiSplit);


if nSmpLabel+nSmpUnlabel ~= nSmp
    error('input error!');
end

if ~isfield(options,'W') 
    options.WOptions.gnd = gnd;
    options.WOptions.semiSplit = semiSplit;
    [W, timeW] = constructW(fea,options.WOptions);
else
    W = options.W;
    timeW = 0;
end

gnd = gnd(semiSplit);

tmp_T = cputime;
   
D = full(sum(W,2));
if isfield(options,'LaplacianNorm') && options.LaplacianNorm
    D=sqrt(1./D);
    D=spdiags(D,0,speye(size(W,1)));
    W=D*W*D;
    L=speye(size(W,1))-W;
else
    L = spdiags(D,0,speye(size(W,1)))-W; 
end

if isfield(options,'LaplacianDegree') 
    L = L^options.LaplacianDegree;
end

elapse.timeW = timeW + cputime - tmp_T;
tmp_T = cputime;

K = constructKernel(fea,[],options);


I=speye(size(K,1));
Ktilde=(I+ReguBeta*K*L)\K;
Ktilde = max(Ktilde,Ktilde');


RegularizerOptions.gnd = gnd;
RegularizerOptions.Kernel = 1;
[eigvector] = KSR_caller(RegularizerOptions, Ktilde(semiSplit,semiSplit));


tempK = I-ReguBeta*L*Ktilde;
eigvector = tempK(:,semiSplit)*eigvector;


if ~isfield(options,'KernelType') || strcmpi(options.KernelType,'Linear') 
    eigvector = fea'*eigvector;
    eigvector = eigvector./repmat(max(1e-10,sum(eigvector.^2,1).^.5),size(eigvector,1),1);
    if nargout >= 3
        Y = fea*eigvector;
    end
else
    tmpNorm = sqrt(sum((eigvector'*K).*eigvector',2));
    eigvector = eigvector./repmat(tmpNorm',size(eigvector,1),1);
    if nargout >= 3
        Y = K*eigvector;
    end
end


elapse.timeMethod = cputime - tmp_T; 
elapse.timeAll =  elapse.timeW + elapse.timeMethod;






