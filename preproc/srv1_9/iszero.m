function [tfFull,tf,n]=iszero(a)
% if vector a is zero vector then tf=true, if a(i)==0, then tfFull(i)=true,
% else =false.
% n: number of zero elements

eps=2^-64;
ab=abs(a);
tfFull=ab<eps;
ma=max(ab);
if ma<eps
   tf=true;
else
    tf=false;
end
n=sum(tfFull);
end