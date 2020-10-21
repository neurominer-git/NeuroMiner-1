function [P, Ppos, nP, Pdesc] = nk_GetModelParams2(analysis, multiflag, CV2ind, curclass )
global MULTI

Pdesc = analysis.Model.ParamDesc{curclass};
if ~exist('curclass','var') || isempty(curclass)
    
    nclass = analysis.SVM.nclass;
    P = cell(nclass,1); Ppos = cell(nclass,1); nP = zeros(nclass,1);
    for curclass = 1 : nclass
        [P{curclass}, Ppos{curclass}, nP(curclass)] = getP(analysis, MULTI, multiflag, CV2ind, curclass);
    end
else
    [P,Ppos, nP] = getP(analysis, MULTI, multiflag, CV2ind, curclass);
end

end

function [P, Ppos, nP] = getP(analysis, MULTI, multiflag, CV2ind, curclass)

if ~isempty(MULTI) && MULTI.flag && multiflag
    % Get optimum multi-group classification parameters
    if iscell(analysis.multi_bestPpos)
        Ppos = analysis.multi_bestPpos{CV2ind};
    else
        Ppos = analysis.multi_bestPpos(CV2ind);
    end
else
    % Get optimum binary classifier / regression parameters 
    if iscell(analysis.bestP{curclass})
        Ppos = analysis.bestPpos{curclass}{CV2ind};
    else
        Ppos = analysis.bestPpos{curclass}(CV2ind);
    end
end

P = analysis.Model.ParamCombs{curclass}(Ppos,:); nP = size(P,1);

end