function [A,Y,numIter,tElapsed,finalResidual]=usr(X,k,option)
% NMF based on multiple update rules: X=AY, s.t. X,A,Y>=0.
% Definition:
%     [A,Y,numIter,tElapsed,finalResidual]=nmfrule(X,k)
%     [A,Y,numIter,tElapsed,finalResidual]=nmfrule(X,k,option)
% X: non-negative matrix, dataset to factorize, each column is a sample, and each row is a feature.
% k: number of clusters.
% option: struct:
% option.distance: distance used in the objective function. It could be
%    'ls': the Euclidean distance (defalut),
%    'kl': KL divergence.
% option.iter: max number of interations. The default is 1000.
% option.dis: boolen scalar, It could be 
%     false: not display information,
%     true: display (default).
% option.residual: the threshold of the fitting residual to terminate. 
%    If the ||X-XfitThis||<=option.residual, then halt. The default is 1e-4.
% option.tof: if ||XfitPrevious-XfitThis||<=option.tof, then halt. The default is 1e-4.
% A: matrix, the basis matrix.
% Y: matrix, the coefficient matrix.
% numIter: scalar, the number of iterations.
% tElapsed: scalar, the computing time used.
% finalResidual: scalar, the fitting residual.
% References:
%  [1]\bibitem{NMF_Lee1999}
%     D.D. Lee and S. Seung,
%     ``Learning the parts of objects by non-negative matrix factorization,''
%     {\it Science},
%     vol. 401, pp. 788-791, 1999.
%  [2]\bibitem{NMF_Lee2001}
%     D.D. Lee and S. Seung,
%     ``Algorithms for non-negative matrix factorization,''
%     {\it Advances in Neural Information Processing Systems}, 
%     vol. 13, pp. 556-562, 2001.
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
% May 01, 2011
%%%%

tStart=tic;
optionDefault.alpha2=1;
optionDefault.alpha1=0.01;
optionDefault.lambda2=0;
optionDefault.lambda1=0.01;
optionDefault.t1=true;
optionDefault.t2=true;
optionDefault.iter=100;
optionDefault.dis=true;
optionDefault.residual=1e-4;
optionDefault.tof=1e-4;
if nargin<3
   option=optionDefault;
else
    option=mergeOption(option,optionDefault);
end

[numFeat,numS]=size(X);

% iter: number of iterations
[r,c]=size(X); % c is # of samples, r is # of features
Y=rand(k,c);
% Y(Y<eps)=0;
Y=max(Y,eps);
A=X/Y;
% A(A<eps)=0;
A=max(A,eps);
XfitPrevious=Inf;
I1=option.alpha2.*eye(k);
I2=option.lambda2.*eye(k);
% E1=option.alpha1.*ones(numFeat,k);
% E2=option.lambda1.*ones(k,numS);
ep=0;
for i=1:option.iter
    H1=Y*Y'+I1;
    if option.t1
        G1=option.alpha1-Y*X';
        A=NNQPActiveSet(H1,G1);
        A=max(A,ep);
    else
        G1=-Y*X';
        A=l1QPActiveSet(H1,G1,option.alpha1);
    end
    A=A';
    %             A(A<eps)=0;
    H2=A'*A+I2;
    if option.t2
        G2=option.lambda1-A'*X;
        Y=NNQPActiveSet(H2,G2);
        Y=max(Y,ep);
    else
        G2=-A'*X;
        Y=l1QPActiveSet(H2,G2,option.lambda1);
    end
    
    if mod(i,10)==0 || i==option.iter
        if option.dis
            disp(['Iterating >>>>>> ', num2str(i),'th']);
        end
        XfitThis=A*Y;
        fitRes=matrixNorm(XfitPrevious-XfitThis);
        XfitPrevious=XfitThis;
        curRes=norm(X-XfitThis,'fro');
        if option.tof>=fitRes || option.residual>=curRes || i==option.iter
            s=sprintf('Mutiple update rules based Sparse NMF on Both Factors successes! \n # of iterations is %0.0d. \n The final residual is %0.4d.',i,curRes);
            disp(s);
            numIter=i;
            finalResidual=curRes;
            break;
        end
    end
end
tElapsed=toc(tStart);
end
