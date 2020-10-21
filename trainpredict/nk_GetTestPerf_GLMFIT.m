function [rs, ds] = nk_GetTestPerf_GLMFIT(~, tXtest, ~, md, ~, ~)
global MODEFL

Z = md.beta(1) + tXtest * md.beta(2:end); 
if strcmp(MODEFL,'classification')
   ds = 1 ./ (1 + exp(-Z)); % Logistic function
else
   ds = Z;
end

ds = ds(:,1);
rs = sign(ds-0.5);



