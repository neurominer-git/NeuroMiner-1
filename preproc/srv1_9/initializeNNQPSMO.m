function [x,s,logV,x1,s1,indx1]=initializeNNQPSMO(H,g)
% initialize the l1QPSMO algorithm
% Yifeng Li
% Mar. 14, 2013

n=numel(g);
x=zeros(n,1);
s=zeros(n,1);

diff=0;
x1=0;
s1=0;
indx1=0;
logV=false;
for i=1:n
    s(i)=g(i);
    if s(i)<-diff
       x1=x(i);
       s1=s(i);
       indx1=i;
       diff=abs(s(i));
       logV=true;
    end
end
end