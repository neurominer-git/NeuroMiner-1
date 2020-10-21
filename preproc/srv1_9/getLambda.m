function Lambda=getLambda(A,G,W)
% solve A_wk*Lambda_wk=G
% Yifeng Li, Feb. 21, 2012

[n2V,nP]=size(W);
% code=2.^(n2V-1:-1:0) * W;
% [codeSorted,csInd]=sort(code);
code=genCode(W);
codeUnik=unique(code);
numUnik=numel(codeUnik);
Lambda=zeros(n2V,nP);
for i=1:numUnik
   indLog=(code==codeUnik(i)); % index of the same working set
   indNum=find(indLog);
   Wi=W(:,indNum(1)); % this unique working set
   if any(~(Wi(1:n2V/2)|Wi(n2V/2+1:end)))
      error('there exist unbounded problem'); 
   end
   Ait=A(Wi,:)';
   Lambdai=Ait\G(:,indNum);
   Lambda(Wi,indNum)=Lambdai;
end
end
