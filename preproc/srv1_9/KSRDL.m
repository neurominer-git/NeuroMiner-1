function [AtA,Y,A,pinvY,numIter,tElapsed,finalResidual]=KSRDL(X,XtX,k,option)
% Kernel dictionary learning, X_{m x n} = A_{m x k} Y_{k x n} given training data X. Each column of X
% is a multivariate sample/signal/instance/observation. A is dictionary/basis matrix/factor loading matrix. Y
% is coefficient matrix / factor score matrix. Note this function only
% returns the AtA instead of A. Therefore if you need A, after calling this function, you can obtain it
% by A=X*pinvY; A=normc(A); if the dictionary prior is uniform. If the
% dictionary prior is Gaussian, you can obtain A via A=X*pinvY.
% XtX: the inner product or kernel matrix of training data X.
% k: the number of dictionary atoms.
% option: struct.
% option.SCMethod: string, specify the spase coding algorithm. It can be
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
% option.dis: logical scalar, if display the optimization process. The default is true.
% option.residual: scalar, regression residual limit to terminate an
% algorithm. The default is 1e-4.
% option.tol: numerical tolerance. The default is 1e-4. 
% AtA: the inner product or kernel matrix of the dictioanary
% Y: the sparse coefficient matrix.
% pinvY: the (regularized) pseudoinverse of Y
% numIter: the number of iterations taken.
% tElapsed: the time taken.
% finalResidual: final regression residual.
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
optionDefault.kernel='linear';
optionDefault.param=[];
optionDefault.dicPrior='uniform';
optionDefault.iter=100;
optionDefault.dis=0;
optionDefault.residual=1e-4;
optionDefault.tof=1e-4;
if nargin<4
    option=optionDefault;
else
    option=mergeOption(option,optionDefault);
end

if strcmp(option.dicPrior,'uniform')
    option.alpha=0;
end
c=size(XtX,1);

if ~isempty(X)
   Xtf=true; 
else
    Xtf=false;
    A=[];
end

% Initialize AtA
% [U,~,~]=svd(X,'econ');
% A=U(:,1:k);
% AtA=computeKernelMatrix(A,A,option);
% AtX=computeKernelMatrix(A,X,option);
% XtX=computeKernelMatrix(X,X,option);
Y=rand(k,c);
% switch option.dicPrior
%     case 'uniform'
%         if c==k
%             AtA=Y'\XtX/Y;
%             AtX=Y'\XtX;
%         else
%             pinvY=pinv(Y);
%             AtA=pinvY'*XtX*pinvY;
%             AtX=pinvY'*XtX;
%         end
%     case 'Gaussian'
%         pinvY=Y'/(Y*Y'+diag( repmat(option.alpha,size(Y,1),1) ) );
%         AtA=pinvY'*XtX*pinvY;
%         AtX=pinvY'*XtX;
% end

prevRes=Inf;
% alternative updating
for i=1:option.iter
    %     % normalize AtA and AtX
    %     [AtA,~,XtX,~,AtX]=normalizeKernelMatrix(AtA,XtX,AtX);
    
   % update kernel matrices
    switch option.dicPrior
        case 'uniform'
            if c==k
                if Xtf
                    A=X/Y;
                end
                AtA=Y'\XtX/Y;
                AtX=Y'\XtX;
            else
                pinvY=pinv(Y);
                if Xtf
                    A=X*pinvY;
                end
                AtA=pinvY'*XtX*pinvY;
                AtX=pinvY'*XtX;
            end
            % normalize AtA and AtX
            if Xtf
               A=normc(A); 
            end
            [AtA,~,XtX,~,AtX]=normalizeKernelMatrix(AtA,XtX,AtX);
        case 'Gaussian'
                pinvY=Y'/(Y*Y'+diag( repmat(option.alpha,size(Y,1),1) ) );
                if Xtf
                    A=X*pinvY;
                end
                AtA=pinvY'*XtX*pinvY;
                AtX=pinvY'*XtX; 
    end
    
    % update Y
    [zeroA,tf0]=iszero(sum(AtA));
    AtA(zeroA,:)=[];
    AtA(:,zeroA)=[];
    AtX(zeroA,:)=[];
    if Xtf
       [zeroA,tf0]=iszero(sum(A,1)); 
       A(:,zeroA)=[];
    end
    if tf0 && k>1
        k=k-1;
    end
    if k<1
       error('k<1, impossible!'); 
    end
    
    Y=KSRSC(AtA,AtX,diag(XtX),option);
%     switch option.SCMethod
%         case 'l1lsIP'    % interior-point method
%             Y=l1QPIPMulti(AtA,AtX,diag(XtX),option);
%         case 'nnlsIP'    % interior-point method
%             Y=NNLSKernelBatchDL(AtA,AtX,diag(XtX),option);
%         case 'l1nnlsIP'  % interior-point method
%             Y=l1NNLSKernelBatchDL(AtA,AtX,diag(XtX),option);
%         case 'l1lsAS'   % active-set method
%             Y=l1QP(AtA,-AtX,option.lambda);
%         case 'nnlsAS'   % active-set method
%             Y=l1NNQPActiveSet(AtA,-AtX,0);
%         case 'l1nnlsAS' % active-set method
%             Y=l1NNQPActiveSet(AtA,-AtX,option.lambda);
%         case 'l1lsPX' % proximal method
%             Y=l1LSProximal(AtA,AtX,diag(XtX),option);
%     end
%     
    
    if mod(i,20)==0 || i==option.iter
        if option.dis
            disp(['Iterating >>>>>> ', num2str(i),'th']);
        end
        switch option.SCMethod
            case {'l1lsIP','l1qpIP','l1lsAS','l1qpAS','l1lsPX','l1qpPX','l1lsSMO','l1qpSMO'}
                curRes=0.5*trace(XtX-Y'*AtX-AtX'*Y+Y'*AtA*Y) + 0.5*option.alpha*trace(AtA) + option.lambda*sum(sum(abs(Y)));
            case {'nnlsIP','nnqpIP','nnlsAS','nnqpAS','nnlsSMO','nnqpSMO'}
                curRes=0.5*trace(XtX-Y'*AtX-AtX'*Y+Y'*AtA*Y) + 0.5*option.alpha*trace(AtA);
        end
        fitRes=prevRes-curRes;
        prevRes=curRes;
        if option.tof>=fitRes || option.residual>=curRes || i==option.iter
            disp(['DL successes!, # of iterations is ',num2str(i),'. The final residual is ',num2str(curRes)]);
            numIter=i;
            finalResidual=curRes;
            break;
        end
    end
end
tElapsed=toc(tStart);
end