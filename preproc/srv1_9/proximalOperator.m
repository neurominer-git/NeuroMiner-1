function x=proximalOperator(xbar,eta)
% solve min_x 1/2||x-xbar||_2^2 + eta||x||_1
% Yifeng Li
% April 11, 2012
% example:
% xbar=[0.2;-0.8;0.3;0.6;-0.1];
% xbar=normc(xbar);
% eta=2^-4;
% x=proximalOperator(xbar,eta)

numx=numel(xbar);
x=xbar;
for i=1:numx
    signxbar=0;
    if xbar(i)>0
        signxbar=1;
    end
    if xbar(i)<0
        signxbar=-1;   
    end
    x(i)=0;
    if abs(xbar(i))>eta
        x(i)=signxbar*(abs(xbar(i))-eta);
    end
end
end