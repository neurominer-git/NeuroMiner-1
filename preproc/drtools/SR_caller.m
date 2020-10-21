function [eigvector, elapse] = SR_caller(options, data)
% SR: Spectral Regression
%
%       [eigvector, elapse] = SR_caller(options, data)
% 
%             Input:
%               data    - data matrix. Each row vector of data is a
%                         sample vector. 
%               options - Struct value in Matlab. The fields in options
%                         that can be set:
%
%                        gnd   -  Colunm vector of the label information
%                                 for each data point. 
%                                 If gnd is provided, SR will give the
%                                 SRDA solution [See Ref 7,8]
%
%                        W     -  Affinity matrix. You can either call
%                                 "constructW" to construct the W, or
%                                 construct it by yourself.
%                                 If gnd is not provided, W is required and
%                                 SR will give the RLPI (RLPP) solution
%                                 [See Ref 1] 
%
%                ReducedDim    -  The number of dimensions. If gnd is
%                                 provided, ReducedDim=c-1 where c is the number
%                                 of classes. Default ReducedDim = 30.
%
%                       Please see SR.m for other options. 
%
%
%             Output:
%               eigvector - Each column is an embedding function, for a new
%                           sample vector (row vector) x,  y = x*eigvector
%                           will be the embedding result of x.
%
%                           If 'Lasso' or 'RidgeLasso' regularization is
%                           used, the output eigvector will be a cell,
%                           please see 'lars.m' for the result format.
%
%               elapse    - Time spent on different steps 
%                           
% 
%===================================================================
%    Examples:
%           
%    (Supervised case with L2-norm (ridge) regularizer, SR-LDA)
%
%       fea = rand(50,70);
%       gnd = [ones(10,1);ones(15,1)*2;ones(10,1)*3;ones(15,1)*4];
%       options = [];
%       options.gnd = gnd;
%       options.ReguAlpha = 0.01;
%       options.ReguType = 'Ridge';
%       [eigvector] = SR_caller(options, fea);
%
%       [nSmp,nFea] = size(fea);
%       if size(eigvector,1) == nFea + 1
%           Y = [fea ones(nSmp,1)]*eigvector;
%       else
%           Y = fea*eigvector;
%       end
%-------------------------------------------------------------------
%    (Unsupervised case with L2-norm (ridge) regularizer, SR-LPP)
%
%       fea = rand(50,70);
%       options = [];
%       options.Metric = 'Euclidean';
%       options.NeighborMode = 'KNN';
%       options.k = 5;
%       options.WeightMode = 'HeatKernel';
%       options.t = 5;
%       W = constructW(fea,options);
%
%       options = [];
%       options.W = W;
%       options.ReguAlpha = 0.01;
%       options.ReguType = 'Ridge';
%       options.ReducedDim = 10;
%       [eigvector] = SR_caller(options, fea);
%
%       [nSmp,nFea] = size(fea);
%       if size(eigvector,1) == nFea + 1
%           Y = [fea ones(nSmp,1)]*eigvector;
%       else
%           Y = fea*eigvector;
%       end
%-------------------------------------------------------------------
%    (Supervised case with L1-norm (lasso) regularizer, SR-SpaseLDA)
%
%       fea = rand(50,70);
%       gnd = [ones(10,1);ones(15,1)*2;ones(10,1)*3;ones(15,1)*4];
%       options = [];
%       options.gnd = gnd;
%       options.ReguAlpha = 0.001;
%       options.ReguType = 'RidgeLasso';
%       options.LassoCardi = [10:5:60];
%       [eigvectorAll] = SR_caller(options, fea);
%       
%       [nSmp,nFea] = size(fea);
%       
%       for i = 1:length(options.LassoCardi)
%           eigvector = eigvectorAll{i};  %projective functions with cardinality options.LassoCardi(i)
%           
%           if size(eigvector,1) == nFea + 1
%               Y = [fea ones(nSmp,1)]*eigvector;
%           else
%               Y = fea*eigvector;
%           end
%       end
%
%-------------------------------------------------------------------
%    (Unsupervised case with L1-norm (lasso) regularizer, SR-SpaseLPP)
%
%       fea = rand(50,70);
%       options = [];
%       options.Metric = 'Euclidean';
%       options.NeighborMode = 'KNN';
%       options.k = 5;
%       options.WeightMode = 'HeatKernel';
%       options.t = 5;
%       W = constructW(fea,options);
%
%       options = [];
%       options.W = W;
%       options.ReguAlpha = 0.001;
%       options.ReguType = 'RidgeLasso';
%       options.LassoCardi = [10:5:60];
%       options.ReducedDim = 10;
%       [eigvector] = SR_caller(options, fea);
%
%       [nSmp,nFea] = size(fea);
%
%       for i = 1:length(options.LassoCardi)
%           eigvector = eigvectorAll{i};  %projective functions with cardinality options.LassoCardi(i)
%           
%           if size(eigvector,1) == nFea + 1
%               Y = [fea ones(nSmp,1)]*eigvector;
%           else
%               Y = fea*eigvector;
%           end
%       end
%
%===================================================================
%
%Reference:
%
%   1. Deng Cai, Xiaofei He, Jiawei Han, "Semi-Supervised Regression using
%   Spectral Techniques", Department of Computer Science
%   Technical Report No. 2749, University of Illinois at Urbana-Champaign
%   (UIUCDCS-R-2007-2749), July 2006.  
%
%   2. Deng Cai, Xiaofei He, Jiawei Han, "Spectral Regression for
%   Dimensionality Reduction", Department of Computer Science
%   Technical Report No. 2856, University of Illinois at Urbana-Champaign
%   (UIUCDCS-R-2007-2856), May 2007.  
%
%   3. Deng Cai, Xiaofei He, Jiawei Han, "SRDA: An Efficient Algorithm for
%   Large Scale Discriminant Analysis", Department of Computer Science
%   Technical Report No. 2857, University of Illinois at Urbana-Champaign
%   (UIUCDCS-R-2007-2857), May 2007.  
%
%   4. Deng Cai, Xiaofei He, and Jiawei Han. "Isometric Projection", Proc.
%   22nd Conference on Artifical Intelligence (AAAI'07), Vancouver, Canada,
%   July 2007.  
%
%   5. Deng Cai, Xiaofei He, Jiawei Han, "Efficient Kernel Discriminant
%   Analysis via Spectral Regression", Department of Computer Science
%   Technical Report No. 2888, University of Illinois at Urbana-Champaign
%   (UIUCDCS-R-2007-2888), August 2007.  
% 
%   6. Deng Cai, Xiaofei He, Jiawei Han, "Spectral Regression: A Unified
%   Subspace Learning Framework for Content-Based Image Retrieval", ACM
%   Multimedia 2007, Augsburg, Germany, Sep. 2007.
%
%   7. Deng Cai, Xiaofei He, Jiawei Han, "Spectral Regression for Efficient
%   Regularized Subspace Learning", IEEE International Conference on
%   Computer Vision (ICCV), Rio de Janeiro, Brazil, Oct. 2007. 
%
%   8. Deng Cai, Xiaofei He, Jiawei Han, "Spectral Regression: A Unified
%   Approach for Sparse Subspace Learning", Proc. 2007 Int. Conf. on Data
%   Mining (ICDM'07), Omaha, NE, Oct. 2007. 
%
%   9. Deng Cai, Xiaofei He, Jiawei Han, "Efficient Kernel Discriminant
%   Analysis via Spectral Regression", Proc. 2007 Int. Conf. on Data
%   Mining (ICDM'07), Omaha, NE, Oct. 2007. 
%
%   10. Deng Cai, Xiaofei He, Wei Vivian Zhang, Jiawei Han, "Regularized
%   Locality Preserving Indexing via Spectral Regression", Proc. 2007 ACM
%   Int. Conf. on Information and Knowledge Management (CIKM'07), Lisboa,
%   Portugal, Nov. 2007.
%
%
%   version 2.0 --Aug/2007 
%   version 1.0 --May/2006 
%
%   Written by Deng Cai (dengcai2 AT cs.uiuc.edu)
%

ReducedDim = 30;
if isfield(options,'ReducedDim')
    ReducedDim = options.ReducedDim;
end

ReducedDim = ReducedDim+1;


[nSmp,nFea] = size(data);

bSup = 0;
if isfield(options,'gnd')
    gnd = options.gnd;
    bSup = 1;
    if length(gnd) ~= nSmp
        error('Label vector wrong!');
    end
else
    if ~isfield(options,'W') | (size(options.W,1) ~= nSmp)
        error('Graph Error!');
    end
    W = options.W;
end

if ReducedDim > nSmp
    ReducedDim = nSmp; 
end


tmp_T = cputime;
if bSup
    Label = unique(gnd);
    nClass = length(Label);

    rand('state',0);
    Y = rand(nClass,nClass);
    Z = zeros(nSmp,nClass);
    for i=1:nClass
        idx = find(gnd==Label(i));
        Z(idx,:) = repmat(Y(i,:),length(idx),1);
    end
    Z(:,1) = ones(nSmp,1);
    [Y,R] = qr(Z,0);
    Y(:,1) = [];
else

    D_mhalf = full(sum(W,2).^-.5);
    if nSmp < 5000
        tmpD_mhalf = repmat(D_mhalf,1,nSmp);
        W = (tmpD_mhalf.*W).*tmpD_mhalf';
        clear tmpD_mhalf;
    else
        [i_idx,j_idx,v_idx] = find(W);
        v1_idx = zeros(size(v_idx));
        for i=1:length(v_idx)
            v1_idx(i) = v_idx(i)*D_mhalf(i_idx(i))*D_mhalf(j_idx(i));
        end
        W = sparse(i_idx,j_idx,v1_idx);
        clear i_idx j_idx v_idx v1_idx
    end
    W = max(W,W');

    dimMatrix = size(W,2);
    if (dimMatrix > 500 & ReducedDim < dimMatrix/10)
        option = struct('disp',0);
        [Y, eigvalue] = eigs(W,ReducedDim,'la',option);
        eigvalue = diag(eigvalue);
    else
        W = full(W);
        [Y, eigvalue] = eig(W);
        eigvalue = diag(eigvalue);

        [junk, index] = sort(-eigvalue);
        eigvalue = eigvalue(index);
        Y = Y(:,index);
        if ReducedDim < length(eigvalue)
            Y = Y(:, 1:ReducedDim);
            eigvalue = eigvalue(1:ReducedDim);
        end
    end
    
    eigIdx = find(abs(eigvalue) < 1e-6);
    eigvalue (eigIdx) = [];
    Y (:,eigIdx) = [];

    nGotDim = length(eigvalue);
    
    
    idx = 1;
    while(abs(eigvalue(idx)-1) < 1e-12)
        idx = idx + 1;
        if idx > nGotDim
            break;
        end
    end
    idx = idx - 1;

    if(idx > 1)  % more than one eigenvector of 1 eigenvalue
        u = zeros(size(Y,1),idx);

        d_m = 1./D_mhalf;
        cc = 1/norm(d_m);
        u(:,1) = cc./D_mhalf;

        bDone = 0;
        for i = 1:idx
            if abs(Y(:,i)' * u(:,1) - 1) < 1e-14
                Y(:,i) = Y(:,1);
                Y(:,1) = u(:,1);
                bDone = 1;
            end
        end

        if ~bDone
            for i = 2:idx
                u(:,i) = Y(:,i);
                for j= 1:i-1
                    u(:,i) = u(:,i) - (u(:,j)' * Y(:,i))*u(:,j);
                end
                u(:,i) = u(:,i)/norm(u(:,i));
            end
            Y(:,1:idx) = u;
        end
    end
    

    if nGotDim < 5000
        Y = repmat(D_mhalf,1,nGotDim).*Y;
    else
        for k = 1:nGotDim
            Y(:,k) = Y(:,k).*D_mhalf;
        end
    end
    
    Y(:,1) = [];
end

elapse.timeResponse = cputime - tmp_T;


if issparse(data)
    options.AppendConstant = 1;
    options.LSQR = 1;
else
    sampleMean = mean(data);
    data = (data - repmat(sampleMean,nSmp,1));
    options.AppendConstant = 0;
    options.LSQR = 0;
end
options.RemoveMean = 0;

[eigvector, elapse.timeReg] = SR(options, Y, data);

elapse.timeAll = elapse.timeResponse + elapse.timeReg;

