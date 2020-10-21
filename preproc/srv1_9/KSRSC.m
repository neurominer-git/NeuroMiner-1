function [X,tElapsed,sparsity]=KSRSC(AtA,AtB,BtB,option)
% kernel sparse coding B=AX given data B and dictionary A. Each column of B
% is a multivariate sample/signal/observation. Each column of A is a
% dictionary atom/basis vector. 
% AtA: matrix, the inner product or kernel matrix of dictionary.
% AtB: matrix, the inner product or kernel matrix between A and B.
% BtB: column vector, BtB(i) is the inner product or kernel value of the
% i-th sample in B.
% option: struct.
% option.SRMethod: string, specify the spase coding algorithm. It can be
% 'l1lsAS': active-set l_1LS sparse coding algorithm
% 'nnlsAS': active-set NNLS sparse coding algorithm (Default)
% 'l1nnlsAS': active-set l_1 NNLS sparse coding algorithm
% 'l1lsIP': interior-point l_1LS sparse coding algorithm
% 'nnlsIP': interior-point NNLS sparse coding algorithm
% 'l1nnlsIP': interior-point l_1 NNLS sparse coding algorithm
% 'l1lsPX': proximal l_1LS sparse coding algorithm.
% option.lambda: non-neagtive scalar, control the parsity of the sparse
% coefficient matrix.
% option.kernel: string, the name of a kernel function, please see function
% computeKernelMatrix function for values of this field. The default is
% 'linear'.
% option.param: column vector, the parameter of the kernel function.
% Default is [].
% option.iter: positive interger, the maximum number of iterations. The
% default is 200.
% option.dis: if display the optimization process. Useless currently
% option.residual: scalar, regression residual limit to terminate an
% algorithm. The default is 1e-4.
% option.tol: numerical tolerance. The default is 1e-4. 
% X: the sparse coefficient matrix.
% tElapsed: time taken.
% sparsity: the sparsity of the coefficient matrix.
%%%%
% Reference:
% [1] Y. Li and A. Ngom, "Sparse representation approaches for the
% classificaiton of high-dimensional biological data," BMC Systems Biology.
% 2013.
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
% Feb. 10, 2012
%%%%

tStart=tic;
optionDefault.lambda=0.1;
optionDefault.SCMethod='nnlsAS';
optionDefault.iter=200;
optionDefault.dis=0;
optionDefault.residual=1e-4;
optionDefault.tof=1e-4;
if nargin<3
    option=optionDefault;
else
    option=mergeOption(option,optionDefault);
end

switch option.SCMethod
    case {'l1lsAS','l1qpAS'}
        X=l1QPActiveSet(AtA,-AtB,option.lambda);
    case {'l1lsIP','l1qpIP'}
        X=l1QPIPMulti(AtA,AtB,BtB,option);
    case {'l1lsPX','l1qpPX'}
        X=l1QPProximalMulti(AtA,AtB,BtB,option);
    case {'l1lsSMO','l1qpSMO'}
        X=l1QPSMOMulti(AtA,-AtB,option.lambda);
    case {'nnlsIP','l1nnlsIP','nnqpIP'}
        X=NNQPIPMulti(AtA,-option.lambda+AtB,BtB,option);
    case {'nnlsAS','l1nnlsAS','nnqpAS'}
        X=NNQPActiveSet(AtA,option.lambda-AtB);
    case {'nnlsSMO','l1nnlsSMO','nnqpSMO'}
        X=NNQPSMOMulti(AtA,option.lambda-AtB);
%     case {'l1nnlsIP','l1nnqpIP'}
%         X=l1NNLSKernelBatchDL(AtA,AtB,BtB,option);
%     case {'l1nnlsAS','l1nnqpAS'}
%         X=l1NNQPActiveSet(AtA,-AtB,option.lambda);
end
sparsity=computeSparsity(X);
tElapsed=toc(tStart);
end