function x=l1QPIP(AtA,Atb,btb,lambda,AtAInv,option)
% Interior-point algorithm to solve the single kernel l1LS/l1QP problem:
% min f(x)=1/2||phi(b)-phi(A)x||_2^2 + lambda||x||_1 where x is of size n
% by 1
% AtA: n by n matrix, the inner product or kernel matrix of A
% Atb: n by 1 matrix, the inner product or kernel matrix between A and b
% btb: scalar, b'*b
% AtAInv: n by n matrix, the inverse of AtA, reserved for future use
% option
% x, n by 1 column vector, the sparse coefficient
% Yifeng Li
% Feb. 07, 2011

if nargin<6
    option=[];
end
optionDefault.NewtonMaxIter=100;
optionDefault.tol=1e-4;
optionDefault.tIni=1/lambda; % the initial value of t
optionDefault.updatetmu=2;
optionDefault.alpha=0.01; % for backtracking line search
optionDefault.beta=0.5; % for backtracking line search
optionDefault.lineSearchMaxIter=5000;
optionDefault.smin=0.5;
optionDefault.pcgTol=1e-3;
optionDefault.maxit=100;
option=mergeOption(option,optionDefault);
nVar=size(AtA,1);
% initialize x and u
x=zeros(nVar,1);%AtAInv*Atb;%zeros(nVar,1);%AtAInv*Atb;
u=ones(nVar,1);
dualGap=Inf;
i=0;
t=option.tIni;
dualVal=-Inf;
while dualGap>=option.tol && i+1<=option.NewtonMaxIter
    i=i+1;
%     fprintf('The %d-th iteration of Newtons method for l1LS ...\n',i);
   % compute gradient and hessian, and then the update deltaxu
   uu_xx=u.*u-x.*x;
   g=[t*(AtA*x-Atb)+2*(x./(u.*u-x.*x));t*lambda*ones(nVar,1)-2*(u./uu_xx)];
   D1=2*((u.*u + x.*x)./(uu_xx.*uu_xx));
   D2=-4*((u.*x)./(uu_xx.*uu_xx));
   H=[t*AtA+diag(D1),diag(D2);diag(D2),diag(D1)];
%    deltaxu=-(H\g);
%    deltaxu=-((H+2^(-32)*eye(2*nVar))\g);
    normg   = norm(g);
    pcgtol  = min(1e-1,option.pcgTol*dualGap/min(1,normg));
    M=[t*diag(diag(AtA))+diag(D1),diag(D2);diag(D2),diag(D1)];
   [deltaxu,flag,relres,iter,resvec]=pcg(H,-g,pcgtol,option.maxit,M);
   % search the step size of Newton's method
    % update (x,u)
   s=1;
   oldf=0.5*t*x'*AtA*x - t*Atb'*x + 0.5*t*btb + t*lambda*sum(u)-sum(log(uu_xx));
   for j=1:option.lineSearchMaxIter
       newx=x+s*deltaxu(1:nVar);
       newu=u+s*deltaxu(nVar+1:2*nVar);
       nunu_nxnx=newu.*newu-newx.*newx;
       if min(nunu_nxnx>0)
           newf1=0.5*t*newx'*AtA*newx - t*Atb'*newx + 0.5*t*btb + t*lambda*sum(newu)-sum(log(nunu_nxnx));
           newf2 = oldf + option.alpha*s*g'*deltaxu;
           if newf1<=newf2
               break;
           end
       else 
%            fprintf('infeasible solution\n');
       end
       s=option.beta*s;
   end
   % update (x,u) use step size s
   x=x+s*deltaxu(1:nVar);
   u=u+s*deltaxu(nVar+1:2*nVar);
   % compute the duality gap
   xtAtAx=x'*AtA*x;
   dualGap=xtAtAx - Atb'*x + lambda*sum(abs(x));
%    dualVal=max(-0.5*xtAtAx+0.5*btb,dualVal);
%    primalVal=0.5*xtAtAx - Atb'*x + 0.5*btb + lambda*sum(abs(x));
%    dualGap=dualGap/dualVal
   % update t
   if s>=option.smin
       t=max(option.updatetmu*min(2*nVar/dualGap,t),t);
   end
%    dualGap=dualGap/dualVal; % relative
end
fprintf('l1LS via Newtons method finished.\n');
end

