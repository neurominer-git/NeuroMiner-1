% =========================================================================
% FORMAT [Y, beta] = nk_PartialCorrelations(G, Y, beta, revertflag)
% =========================================================================
% Remove nuisance effects G from Y (optionally, using a predefined beta)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Nikolaos Koutsouleris, 07 / 2011
function [Y, beta] = nk_Regress(G, Y, beta, revertflag)

%fprintf(1,'\n* Remove covariable effects.')
if ~exist('beta','var') || isempty(beta), 
    beta = zeros(size(Y,2), size(G,2));
    T = zeros(size(Y));
    for i=1:size(G,2)
        for j=1:size(Y,2)
            [tbeta, dum, T(:,j) ] = regress(Y(:,j), [ones(size(G,1),1) G(:,i)]);
            beta(j,i) = tbeta(2);
        end
    end
end;

return
