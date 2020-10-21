function DScore = simbamain(Y, labels, cpus, SortInd)

global FEATSEL
gpu = FEATSEL.simba.gpu;
extra_param = FEATSEL.simba.extra_param;

if nargin < 5, if ~isempty(SortInd), Y = resamp(Y, SortInd, cpus); end; end
if FEATSEL.salthreshmode == 1, extra_param.salCI = FEATSEL.salCI; end

if isfield(extra_param,'beta') && strcmp(extra_param.beta,'auto'), 
    extra_param.beta = suggestBeta(Y, labels);
    fprintf(' (beta=%g)',extra_param.beta)
end
DScore = nk_SimbaMain(Y, labels, extra_param, gpu);