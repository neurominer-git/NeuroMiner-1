function cost=IMRelief_cost(alpha, Z, Weight, descent, Targets, lambda)

%Function ComputeAlpha: compute alpha using linear search

Weight = Weight-alpha*descent;
a = ((Weight(:).^2)')*Z; %margin
Result = 1./(1+exp(-a)); %probability of being class 1

index = Result==1;
Result(index) = 1-10^(-10);
index = Result==0;
Result(index) = 10^(-10);
 
if sum(Result==1)>=1 || sum(Result==0)>=1
    cost = +inf;
else
    cost = -sum(Targets.*log(Result(:))+(1-Targets).*log(1-Result(:)))+lambda*((Weight(:))'*Weight(:));
end

return
