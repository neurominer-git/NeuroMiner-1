function [Y,numIter,tElapsed]=kernelseminmfruletest(X,outTrain)
% UR-NNLS: X=AY, s.t. Y>0;
% X, matrix, each column is a new sample, each row is a feature
% outTrain: struct, with two fiels:
% outTrain.factors{1} is the dictonary A.
% outTrain.option: struct, the option to setup, with fields:
% option.iter, max number of interations
% option.kernel, string, kernel function
% option.param, scalar or column vector, the parameter of kernel
% option.dis, boolen scalar, false: (not display information) or true
% (display).
% option.residual, if the ||X-XfitThis||<=option.residual, halt.
% option.tof, if ||XfitPrevious-XfitThis||<=option.tof, halt.
% Yifeng Li, May 01, 2011
% yifeng.li.cn@gmail.com

tStart=tic;
optionDefault.iter=1000;
optionDefault.kernel='linear';
optionDefault.param=[];
optionDefault.dis=1;
optionDefault.residual=1e-4;
optionDefault.tof=1e-4;

A=outTrain.factors{1};
if isfield(outTrain,'option')
    option=outTrain.option;
else
    option=[];
end
option=mergeOption(option,optionDefault);

% compute kernels
Apn=computeKernelMatrix(A,X,option);
Bpn=computeKernelMatrix(A,A,option);
KX=computeKernelMatrix(X,X,option);
Y=pinv(Bpn)*Apn;

% Y=A\X;
Y(Y<0)=0;
%Apn=A'*X;
Ap=(abs(Apn)+Apn)./2;
An=(abs(Apn)-Apn)./2;
%Bpn=A'*A;
Bp=(abs(Bpn)+Bpn)./2;
Bn=(abs(Bpn)-Bpn)./2;
prevRes=Inf;
for i=1:option.iter
    Y=Y.*sqrt((Ap + Bn*Y)./(An + Bp*Y));
    if mod(i,100)==0 || i==option.iter
        if option.dis
            disp(['Iterating >>>>>> ', num2str(i),'th']);
        end
        YtApn=Y'*Apn;
        curRes=trace(KX-YtApn-YtApn'+Y'*Bpn*Y); % use abs because the trace maybe negative due to numerical reasons
        fitRes=prevRes-curRes;%abs(prevRes-curRes);
        prevRes=curRes;
%         XfitThis=A*Y;
%         fitRes=matrixNorm(XfitPrevious-XfitThis);
%         XfitPrevious=XfitThis;
%         curRes=matrixNorm(X-XfitThis);
        if option.tof>=fitRes || option.residual>=curRes || i==option.iter
            disp(['UR-NNLS successes!, # of iterations is ',num2str(i),'. The final residual is ',num2str(curRes)]);
            numIter=i;
            finalResidual=curRes;
            break;
        end
    end
end
tElapsed=toc(tStart);
end
