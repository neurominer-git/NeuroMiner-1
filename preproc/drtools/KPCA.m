function [eigvector, eigvalue, elapse, Y] = KPCA(options, data, ReducedDim)
%KPCA	Kernel Principal Component Analysis
%
%	Usage:
%       [eigvector, eigvalue] = KPCA(options, data, ReducedDim)
%       [eigvector, eigvalue] = KPCA(options, data)
% 
%             Input:
%             options   - Struct value in Matlab. The fields in options
%                         that can be set:
%                      Kernel  -  1: data is actually the kernel matrix. 
%                                 0: ordinary data matrix. 
%                                   Default: 0 
%                         
%                       Please see constructKernel.m for other Kernel options. 
%
%               data    - 
%                      if options.Kernel = 0
%                           Data matrix. Each row vector of fea is a data
%                           point. 
%                      if options.Kernel = 1
%                           Kernel matrix. 
%
%          ReducedDim   - The dimensionality of the reduced subspace. If 0,
%                         all the dimensions will be kept. 
%                         Default is 0. 
%
%             Output:
%               eigvector - Each column is an embedding function, for a new
%                           data point (row vector) x,  y = K(x,:)*eigvector
%                           will be the embedding result of x.
%                           K(x,:) = [K(x1,x),K(x2,x),...K(xm,x)]
%               eigvalue  - The sorted eigvalue of PCA eigen-problem. 
%
%	Examples:
%           options.KernelType = 'Gaussian';
%           options.t = 1;
% 			fea = rand(7,10);
% 			[eigvector,eigvalue] = KPCA(options,fea,4);
%           feaTest = rand(3,10);
%           Ktest = constructKernel(feaTest,fea,options)
%           Y = Ktest*eigvector;
% 
%Reference:
%
%   Bernhard Schölkopf, Alexander Smola, Klaus-Robert Müller, “Nonlinear
%   Component Analysis as a Kernel Eigenvalue Problem", Neural Computation,
%   10:1299-1319, 1998.
%
% 
%   version 1.0 --April/2005 
%
%   Written by Deng Cai (dengcai2 AT cs.uiuc.edu)
%                                                   


if (~exist('ReducedDim','var'))
   ReducedDim = 0;
end

if isfield(options,'Kernel') & options.Kernel
    K = data;
    clear data;
    K = max(K,K');
    elapse.timeK = 0;
else
    [K, elapse.timeK] = constructKernel(data,[],options);
end
    
nSmp = size(K,1);
if (ReducedDim > nSmp) | (ReducedDim <=0)
    ReducedDim = nSmp;
end

tmp_T = cputime;

K_old = K;

sumK = sum(K,2);
H = repmat(sumK./nSmp,1,nSmp);
K = K - H - H' + sum(sumK)/(nSmp^2);
K = max(K,K');
clear H;

if nSmp > 1000 & ReducedDim < nSmp/10  
    % using eigs to speed up!
    option = struct('disp',0);
    [eigvector, eigvalue] = eigs(K,ReducedDim,'la',option);
    eigvalue = diag(eigvalue);
else
    [eigvector, eigvalue] = eig(K);
    eigvalue = diag(eigvalue);
    
    [junk, index] = sort(-eigvalue);
    eigvalue = eigvalue(index);
    eigvector = eigvector(:,index);
end


if ReducedDim < length(eigvalue)
    eigvalue = eigvalue(1:ReducedDim);
    eigvector = eigvector(:, 1:ReducedDim);
end

maxEigValue = max(abs(eigvalue));
eigIdx = find(abs(eigvalue)/maxEigValue < 1e-6);
eigvalue (eigIdx) = [];
eigvector (:,eigIdx) = [];

for i=1:length(eigvalue) % normalizing eigenvector
    eigvector(:,i)=eigvector(:,i)/sqrt(eigvalue(i));
end;


elapse.timeMethod = cputime - tmp_T; 
elapse.timeAll = elapse.timeK + elapse.timeMethod;

if nargout >= 4
    Y = K_old*eigvector;
end
