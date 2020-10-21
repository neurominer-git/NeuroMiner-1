function [varind, varstr] = nk_SelectCovariateIndex(dat, varind, askfl)

if ~exist('varind','var') || isempty(varind), varind = 1; end
if ~exist('askfl','var') || isempty(askfl), askfl = 1; end

nD = size(dat.covars,2);
if askfl && nD == 1, 
    varind = 1;
else
    fprintf('\n******************************')
    fprintf('\n*****    COVARIATE(S)    *****')  
    fprintf('\n******************************')

    for i=1: nD
        fprintf('\n(%g)\t%s',i,dat.covnames{i});
    end

    if askfl  && nargout > 0, varind = nk_input('Select covariate(s)',0,'i', varind); end
end
varstr = ['_var' num2str(varind)];

return