function [besta,bestalpha,bestb,bestMin]=learnCurve(es,ns,option)
% learn the learning curve
% es: vector, the error rates
% ns: vector, the correponding numbers of trining samples
% option: struct, reserved, not used now
% Yifeng Li
% Oct 29, 2012
% example:
% es=[0.4;0.28;0.17;0.04;0.01;0.005;0.0005;0.0005];
% ns=[10; 20;  30; 40;50;60;80;100];
% [besta,bestalpha,bestb,bestMin]=learnCurve(es,ns)

if nargin<3
   option=[]; 
end
optionDefault.someField=[];
option=mergeOption(option,optionDefault);


f=@(x) invPowerLaw(x,ns,es);
x0=[0;1;0.01];
lb=[-100;0;0];
ub=[100;10;1];
[x,resnorm] = lsqnonlin(f,x0,lb,ub);
besta=x(1);
bestalpha=x(2);
bestb=x(3);
bestMin=resnorm;

% 
% bs=(0.2:-option.bstep:0)';
% numb=numel(bs);
% minRs=nan(numb,1);
% as=nan(numb,1);
% alphas=nan(numb,1);
% bestb=0;
% bestMin=inf;
% besta=0;
% bestalpha=0;
% for i=1:numb
% %     d=log(b-es);
% %     C=[-ones(numel(ns),1),ns];
% %     lb=[-100;0];
% %     ub=[];
% %     Aeq=[];
% %     beq=[];
% %     A=[];
% %     bb=[];
% %     [x,resnorm] = lsqlin(C,d,A,bb,Aeq,beq,lb,ub);
%     as(i)=exp(x(1));
%     alphas(i)=x(2);
%     minRs(i)=resnorm;
% end
% [bestMin,ind]=min(minRs);
% besta=as(ind);
% bestalpha=alphas(ind);
% bestb=bs(i);
end
