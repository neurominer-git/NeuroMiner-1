function K=dynamicSystems_kernel(D1,D2,param)
% Dynamical system kernel function
% D1,D2, matrice, vectorized samples, each column is a sample
% Reference:
% K.M. Borgwardt, Class prediction from time series gene expression
% profiles using dynamical systems kernels. Pracific Symposium on
% Biocomputing, vol. 11, pp. 547-558, 2006.


numR=param(1);
numC=param(2);
rank=param(3);
lambda=param(4);
D1=D1';
D2=D2';
numP1=size(D1,1); % number of pages
numP2=size(D2,1); % number of pages
D1=unmatrizicing(D1,3,[numR,numC,numP1]);
D2=unmatrizicing(D2,3,[numR,numC,numP2]);
K=nan(numP1,numP2);
for p1=1:numP1
    % estimate the parameters (P,Q,R,SS) of sample p1
   [P1,Q1,R1,S1,X1]=estimatePQRS(D1(:,:,p1),rank);   
    for p2=1:numP2
         % estimate the parameters (P,Q,R,S) of sample p2
       [P2,Q2,R2,S2,X2]=estimatePQRS(D2(:,:,p2),rank);
         % calculate M1 and M2
         M1=dlyap(exp(-lambda)*Q1',Q2',exp(-lambda)*Q1'*P1'*P2*Q2);
         M2=dlyap(exp(-lambda)*Q1',Q2',P1'*P2);
         %calculate the similarity of samples p1 and p2
        K(p1,p2)=X1(:,1)'*M1*X2(:,1) + (1/(exp(lambda) -1))*(trace(S1*M2)+trace(R1));
%         if p2<numP1
%            K(p2,p1)=K(p1,p2); 
%         end
    end
end
end