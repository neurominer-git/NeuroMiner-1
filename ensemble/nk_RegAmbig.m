function param = nk_RegAmbig(sE, tE, ambigtype)
% =========================================================================
% FORMAT function param = nk_RegAmbig(E)
% =========================================================================
%
% this function measures the average divergence of single ensemble component 
% predictions from the ensemble prediction
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Nikolaos Koutsouleris, 10/2011

param = 0; 
if isempty(sE), return; end
if ~exist('tE','var'), 
    tE=[];
    if ~exist('ambigtype','var'), ambigtype='regambig'; end
else
    if ~exist('ambigtype','var'), ambigtype='pearson'; end 
end

switch ambigtype

    case 'regambig'
        E = [sE tE];
        mE = median(E,2);
        param = mean(sum((bsxfun(@minus, E, mE ).^2)));

    otherwise
        param = 1 - abs(nk_CorrMat( median(sE,2),tE,ambigtype));

end