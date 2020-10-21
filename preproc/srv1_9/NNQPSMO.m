function [x,s,numIter]=NNQPSMO(H,g)
% The SMO algorithm for l1QP problem: 0.5*x'Hx + g'x s.t. x>=0
% Yifeng Li
% Mar. 17, 2013

tol=1e-3;
n=numel(g);
numIter=0;
[x,s,logV,x1,s1,indx1]=initializeNNQPSMO(H,g);
if ~logV
    return;
end

% [logV,s1,x1,indx1]=findViolateXl1QPSMO(s,x,lambda);
while logV
    x1Old=x1;
    s1Old=s1;
    %     xOld=x;
    sOld=s;
    h11=H(indx1,indx1);
    xnp=(-s1Old+h11*x1Old)/h11;
    if xnp<0
        x1=0;
    else
        x1=xnp;
    end
    % update x,s
    x(indx1)=x1;
    xd=x1-x1Old;
%     xmax=max(x);
%     xtol=abs(xmax*tol);
    xtol=1e-10;
    diff=1e-10;
    x11=0;
    s11=0;
    indx11=0;
    logV=false;
    for i=1:n
        s(i)=H(i,indx1)*xd + sOld(i);
        % find the best variable which violate the optimality condition
        if -xtol<=x(i) && x(i)<=xtol
            if s(i)<-diff
                x11=x(i);
                s11=s(i);
                indx11=i;
                diff=abs(s(i));
                logV=true;
            end
        else % if x(i)>xtol 
            if abs(s(i))> diff
                x11=x(i);
                s11=s(i);
                indx11=i;
                diff=abs(s(i));
                logV=true;
            end
        end
    end
    x1=x11;
    s1=s11;
    indx1=indx11;
    %     [logV,s1,x1,indx1]=findViolateXl1QPSMO(s,x,lambda);
    numIter=numIter+1;
end
end