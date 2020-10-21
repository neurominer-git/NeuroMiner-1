function [xhat,numIter,finalResidual]=l1QPProximal(AtA,Atb,btb,option)
% Proximal method for single l1ls/l1QP problem
% Yifeng Li
% April 11, 2012

k=size(AtA,1);
if k==0
   error('k==0!'); 
end
optionDefault.lambda=2^-4;
optionDefault.L=k; % set k if sparse coding, set larger such as 1e+4*k if dictionary learning
optionDefault.iter=10*k;
optionDefault.residual=1e-7;
optionDefault.tof=1e-7;
optionDefault.dis=false;
if nargin<4
   option=[]; 
end
option=mergeOption(option,optionDefault);

xhat=zeros(k,1);
prevRes=Inf;
for i=1:option.iter
    xbar=xhat-(AtA*xhat-Atb)/option.L;
    eta=option.lambda/option.L;
    xhat=proximalOperator(xbar,eta);
    if mod(i,10)==0 || i==option.iter
        if option.dis
%             disp(['Iterating >>>>>> ', num2str(i),'th']);
        end
        curRes=0.5*(xhat'*AtA*xhat) -Atb'*xhat + 0.5*btb + option.lambda*(sum(abs(xhat))); 
        %curRes;
        fitRes=prevRes-curRes;
        prevRes=curRes;
        if option.tof>=fitRes || option.residual>=curRes || i==option.iter
            if option.dis
                disp(['Proximal method successes!, # of iterations is ',num2str(i),'. The final residual is ',num2str(curRes)]);
            end
            numIter=i;
            finalResidual=curRes;
            break;
        end
    end 
end
end