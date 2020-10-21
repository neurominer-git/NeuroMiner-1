function P = nk_CalibrateProbabilities(P, cutoff)

if ~exist('cutoff','var') || isempty(cutoff)
    P = P - 0.05; 
else
    P = P - (cutoff+realmin);
end

end