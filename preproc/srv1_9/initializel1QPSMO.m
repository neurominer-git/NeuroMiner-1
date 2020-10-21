function [x,s,logV,x1,s1,indx1]=initializel1QPSMO(H,g,lambda)
% initialize the l1QPSMO algorithm
% Yifeng Li
% Mar. 14, 2013

tol=1e-3;
n=numel(g);
% x=2*lambda*rand(n,1)-lambda;
x=zeros(n,1);
% s=H*x+g;
s=zeros(n,1);
dtol=lambda*tol;

diff=0;
x1=0;
s1=0;
indx1=0;
logV=false;
for i=1:n
    s(i)=g(i);
    diffi=abs(s(i))-lambda;
    if diffi>diff && diffi > dtol
       x1=x(i);
       s1=s(i);
       indx1=i;
       diff=diffi;
       logV=true;
    end
end
end