function [eigvector, eigvalue, elapse] = LGE(W, D, options, data)
% LGE: Linear Graph Embedding
%
%       [eigvector, eigvalue] = LGE(W, D, options, data)
% 
%             Input:
%               data       - data matrix. Each row vector of data is a
%                         sample vector. 
%               W       - Affinity graph matrix. 
%               D       - Constraint graph matrix. 
%                         LGE solves the optimization problem of 
%                         a* = argmax (a'data'WXa)/(a'data'DXa) 
%                         Default: D = I 
%
%               options - Struct value in Matlab. The fields in options
%                         that can be set:
%
%                     ReducedDim   -  The dimensionality of the reduced
%                                     subspace. If 0, all the dimensions
%                                     will be kept. Default is 30. 
%
%                            Regu  -  1: regularized solution, 
%                                        a* = argmax (a'X'WXa)/(a'X'DXa+ReguAlpha*I) 
%                                     0: solve the sinularity problem by SVD (PCA) 
%                                     Default: 0 
%
%                        ReguAlpha -  The regularization parameter. Valid
%                                     when Regu==1. Default value is 0.1. 
%
%                            ReguType  -  'Ridge': Tikhonov regularization
%                                         'Custom': User provided
%                                                   regularization matrix
%                                          Default: 'Ridge' 
%                        regularizerR  -   (nFea x nFea) regularization
%                                          matrix which should be provided
%                                          if ReguType is 'Custom'. nFea is
%                                          the feature number of data
%                                          matrix
%
%                            PCARatio     -  The percentage of principal
%                                            component kept in the PCA
%                                            step. The percentage is
%                                            calculated based on the
%                                            eigenvalue. Default is 1
%                                            (100%, all the non-zero
%                                            eigenvalues will be kept.
%                                            If PCARatio > 1, the PCA step
%                                            will keep exactly PCARatio principle
%                                            components (does not exceed the
%                                            exact number of non-zero components).  
%
%
%             Output:
%               eigvector - Each column is an embedding function, for a new
%                           sample vector (row vector) x,  y = x*eigvector
%                           will be the embedding result of x.
%               eigvalue  - The sorted eigvalue of the eigen-problem.
%               elapse    - Time spent on different steps 
%                           
% 
%
%    Examples:
%
% See also LPP, NPE, IsoProjection, LSDA.
%
%
%Reference:
%
%   Deng Cai, Xiaofei He, Yuxiao Hu, Jiawei Han, and Thomas Huang, 
%   "Learning a Spatially Smooth Subspace for Face Recognition", CVPR'2007
% 
%   version 2.1 --June/2007 
%   version 2.0 --May/2007 
%   version 1.0 --Sep/2006 
%
%   Written by Deng Cai (dengcai2 AT cs.uiuc.edu)
%

if ~exist('data','var')
    global data;
end


if (~exist('options','var'))
   options = [];
end

if isfield(options,'ReducedDim')
    Dim = options.ReducedDim;
else
    Dim = 30;
end

if ~isfield(options,'Regu') | ~options.Regu
    bPCA = 1;
    if ~isfield(options,'PCARatio')
        options.PCARatio = 1;
    end
else
    bPCA = 0;
    if ~isfield(options,'ReguType')
        options.ReguType = 'Ridge';
    end
    if ~isfield(options,'ReguAlpha')
        options.ReguAlpha = 0.1;
    end
end

bD = 1;
if ~exist('D','var') | isempty(D)
    bD = 0;
end


[nSmp,nFea] = size(data);
if size(W,1) ~= nSmp
    error('W and data mismatch!');
end
if bD & (size(D,1) ~= nSmp)
    error('D and data mismatch!');
end

tmp_T = cputime;

bChol = 0;
if bPCA & (nSmp > nFea) & (options.PCARatio >= 1)
    if bD
        DPrime = data'*D*data;
    else
        DPrime = data'*data;
    end
    if issparse(DPrime)
        DPrime = full(DPrime);
    end
    DPrime = max(DPrime,DPrime');
    [R,p] = chol(DPrime);
    
    if p == 0
        bPCA = 0;
        bChol = 1;
    end
end

%======================================
% SVD
%======================================

if bPCA    
    if nSmp > nFea
        ddata = data'*data;
        if issparse(ddata)
            ddata = full(ddata);
        end
        ddata = max(ddata,ddata');

        [eigvector_PCA, eigvalue_PCA] = eig(ddata);
        eigvalue_PCA = diag(eigvalue_PCA);
        clear ddata;

        maxEigValue = max(abs(eigvalue_PCA));
        eigIdx = find(eigvalue_PCA/maxEigValue < 1e-12);
        eigvalue_PCA(eigIdx) = [];
        eigvector_PCA(:,eigIdx) = [];

        [junk, index] = sort(-eigvalue_PCA);
        eigvalue_PCA = eigvalue_PCA(index);
        eigvector_PCA = eigvector_PCA(:, index);
        
        %=======================================
        if options.PCARatio > 1
            idx = options.PCARatio;
            if idx < length(eigvalue_PCA)
                eigvalue_PCA = eigvalue_PCA(1:idx);
                eigvector_PCA = eigvector_PCA(:,1:idx);
            end
        elseif options.PCARatio < 1
            sumEig = sum(eigvalue_PCA);
            sumEig = sumEig*options.PCARatio;
            sumNow = 0;
            for idx = 1:length(eigvalue_PCA)
                sumNow = sumNow + eigvalue_PCA(idx);
                if sumNow >= sumEig
                    break;
                end
            end
            eigvalue_PCA = eigvalue_PCA(1:idx);
            eigvector_PCA = eigvector_PCA(:,1:idx);
        end
        %=======================================
        
        if bD
            data = data*eigvector_PCA;
        else
            eigvalue_PCA = eigvalue_PCA.^-.5;
            data = (data*eigvector_PCA).*repmat(eigvalue_PCA',nSmp,1);
        end
    else
        ddata = data*data';
        if issparse(ddata)
            ddata = full(ddata);
        end
        ddata = max(ddata,ddata');

        [eigvector, eigvalue_PCA] = eig(ddata);
        eigvalue_PCA = diag(eigvalue_PCA);
        clear ddata;

        maxEigValue = max(eigvalue_PCA);
        eigIdx = find(eigvalue_PCA/maxEigValue < 1e-12);
        eigvalue_PCA(eigIdx) = [];
        eigvector(:,eigIdx) = [];

        [junk, index] = sort(-eigvalue_PCA);
        eigvalue_PCA = eigvalue_PCA(index);
        eigvector = eigvector(:, index);
        
        %=======================================
        if options.PCARatio > 1
            idx = options.PCARatio;
            if idx < length(eigvalue_PCA)
                eigvalue_PCA = eigvalue_PCA(1:idx);
                eigvector = eigvector(:,1:idx);
            end
        elseif options.PCARatio < 1
            sumEig = sum(eigvalue_PCA);
            sumEig = sumEig*options.PCARatio;
            sumNow = 0;
            for idx = 1:length(eigvalue_PCA)
                sumNow = sumNow + eigvalue_PCA(idx);
                if sumNow >= sumEig
                    break;
                end
            end
            eigvalue_PCA = eigvalue_PCA(1:idx);
            eigvector = eigvector(:,1:idx);
        end
        %=======================================
        
        eigvalue_PCA = eigvalue_PCA.^.5;
        eigvalue_PCAMinus = eigvalue_PCA.^-1;

        eigvector_PCA = (data'*eigvector).*repmat(eigvalue_PCAMinus',nFea,1);

        if bD
            data = eigvector.*repmat(eigvalue_PCA',nSmp,1);
        else
            data = eigvector;
        end
        
        eigvalue_PCA = eigvalue_PCAMinus;
        

        clear eigvector;
    end
    
    if bD
        DPrime = data'*D*data;
        DPrime = max(DPrime,DPrime');
    end
else
    if ~bChol
        if bD
            DPrime = data'*D*data;
        else
            DPrime = data'*data;
        end

        switch lower(options.ReguType)
            case {lower('Ridge')}
                for i=1:size(DPrime,1)
                    DPrime(i,i) = DPrime(i,i) + options.ReguAlpha;
                end
            case {lower('Tensor')}
                DPrime = DPrime + options.ReguAlpha*options.regularizerR;
            case {lower('Custom')}
                DPrime = DPrime + options.ReguAlpha*options.regularizerR;
            otherwise
                error('ReguType does not exist!');
        end

        DPrime = max(DPrime,DPrime');
    end
end

WPrime = data'*W*data;
WPrime = max(WPrime,WPrime');


elapse.timePCA = cputime - tmp_T;

tmp_T = cputime;


%======================================
% Generalized Eigen
%======================================

dimMatrix = size(WPrime,2);

if Dim > dimMatrix
    Dim = dimMatrix; 
end


if isfield(options,'bEigs')
    if options.bEigs
        bEigs = 1;
    else
        bEigs = 0;
    end
else
    if (dimMatrix > 1000 & Dim < dimMatrix/10) | (dimMatrix > 500 & Dim < dimMatrix/20) | (dimMatrix > 250 & Dim < dimMatrix/30) 
        bEigs = 1;
    else
        bEigs = 0;
    end
end


if bEigs
    %disp('use eigs to speed up!');
    option = struct('disp',0);
    if bPCA & ~bD
        [eigvector, eigvalue] = eigs(WPrime,Dim,'la',option);
    else
        if bChol
            option.cholB = 1;
            [eigvector, eigvalue] = eigs(WPrime,R,Dim,'la',option);
        else
            [eigvector, eigvalue] = eigs(WPrime,DPrime,Dim,'la',option);
        end
    end
    eigvalue = diag(eigvalue);
else
    if bPCA & ~bD 
        [eigvector, eigvalue] = eig(WPrime);
    else
        [eigvector, eigvalue] = eig(WPrime,DPrime);
    end
    eigvalue = diag(eigvalue);
    
    [junk, index] = sort(-eigvalue);
    eigvalue = eigvalue(index);
    eigvector = eigvector(:,index);

    if Dim < size(eigvector,2)
        eigvector = eigvector(:, 1:Dim);
        eigvalue = eigvalue(1:Dim);
    end
end



if bPCA
    if bD
        eigvector = eigvector_PCA*eigvector;
    else
        eigvector = eigvector_PCA*(repmat(eigvalue_PCA,1,length(eigvalue)).*eigvector);
    end
end

for i = 1:size(eigvector,2)
    eigvector(:,i) = eigvector(:,i)./norm(eigvector(:,i));
end

    
elapse.timeMethod = cputime - tmp_T; 
elapse.timeAll = elapse.timePCA + elapse.timeMethod;


