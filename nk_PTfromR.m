function [p, t] = nk_PTfromR(r, v, fl) 

r(isnan(r))=0;
t = r .* sqrt((v-2) ./ (1 - r.^2));
if ~exist('v','var') || isempty(v)
    error('Degrees of freedom missing!')
end
if ~exist('v','var') || isempty(v)
    error('One or two-sided p values have to be specified!')
end
x=1-r.^2; p=betainc(x,0.5*(v-2),0.5);
if fl==1, p = 1-(1-p)/2; end
    