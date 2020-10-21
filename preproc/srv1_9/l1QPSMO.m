function [x,numIter]=l1QPSMO(H,g,lambda)
% The SMO algorithm for l1QP problem: 0.5*x'Hx + g'x + lambda|X|_1
% Yifeng Li
% Mar. 14, 2013

% tol=1e-3;
% dtol=lambda*tol;
xtol=1e-10;
dtol=1e-10;
n=numel(g);
numIter=0;
[x,s,logV,x1,s1,indx1]=initializel1QPSMO(H,g,lambda);
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
    %     xpos=(-s1Old+h11*x1Old-lambda)/h11;
    %     if xpos<0
    %        xpos=0;
    %     end
    %     fxpos=0.5*h11*xpos*xpos + (s1Old-h11*x1Old+lambda)*xpos;
    %     xneg=(-s1Old+h11*x1Old+lambda)/h11;
    %     if xneg>0
    %        xneg=0;
    %     end
    %     fxneg=0.5*h11*xneg*xneg + (s1Old-h11*x1Old-lambda)*xneg;
    %     if fxpos<fxneg
    %        x1=xpos;
    %     else
    %         x1=xneg;
    %     end
    b1=s1Old-h11*x1Old;
    if b1>=0
        sign=1;
    else
        sign=-1;
    end
    if abs(b1)>=lambda
        x1= sign*(lambda-abs(b1))/h11;
    else
        x1=0;
    end
    % update x,s
    x(indx1)=x1;
    xd=x1-x1Old;
    %     xmax=max(x);
    %     xtol=abs(xmax*tol);
    diff=dtol;
    x11=0;
    s11=0;
    indx11=0;
    logV=false;
    for i=1:n
        s(i)=H(i,indx1)*xd + sOld(i);
        % find the best variable which violate the optimality condition
        if x(i)>xtol  %&& abs(s(i)+lambda)> dtol
            if abs(s(i)+lambda)> diff
                diff=abs(s(i)+lambda);
                x11=x(i);
                s11=s(i);
                indx11=i;
                logV=true;
            end
        else
            if x(i)<-xtol %&& abs(s(i)-lambda)> dtol
                if abs(s(i)-lambda)>diff
                    diff=abs(s(i)-lambda);
                    x11=x(i);
                    s11=s(i);
                    indx11=i;
                    logV=true;
                end
            else
                if x(i)<=xtol && x(i)>=-xtol  %&& (abs(s1)-lambda)> dtol
                    if (abs(s(i))-lambda)>diff
                        diff=abs(s(i))-lambda;
                        x11=x(i);
                        s11=s(i);
                        indx11=i;
                        logV=true;
                    end
                end
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