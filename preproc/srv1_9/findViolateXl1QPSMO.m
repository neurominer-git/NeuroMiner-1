function [logV,s1,x1,indx1]=findViolateXl1QPSMO(s,x,lambda)
% find a variable which violating the optimality condition
% Yifeng Li
% Mar. 14, 2013

tol=1e-3;
xmax=max(x);
xtol=abs(xmax*tol);
dtol=lambda*tol;
n=numel(x);
i=randi(n);
logV=false;
for j=1:n
   x1=x(i);
   s1=s(i);
   indx1=i;
   if x1>xtol && abs(s1+lambda)> dtol
       logV=true;
       return;
   end
  if x1<-xtol && abs(s1-lambda)> dtol
       logV=true;
       return;
  end 
   if x1<=xtol && x1>=-xtol && (abs(s1)-lambda)> dtol
       logV=true;
       return;
   end
   i=i+1;
   if i>n
      i=1; 
   end
end
end
