function [P, I] = nk_ComputeEnsembleCasePropagationProbability(C, nNodes)

nCases = size(C,1);
P = zeros(nCases,nNodes);
I = zeros(nCases,1);
for i=1:nCases
    ni = numel(C{i});
    for j=1:nNodes
        P(i,j) = sum(C{i}==j) * 100 / ni;
    end
    [~,I(i)] = max(P(i,:),[],2);
end