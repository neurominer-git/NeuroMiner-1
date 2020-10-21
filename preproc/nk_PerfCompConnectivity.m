function [ C, params ] = nk_PerfCompConnectivity( Y, params )

% Step 1: Get Seed point for connectivity analysis if not specified
if ~isfield(params,'seeds')
    if isfield(params,'seedparams') 
        fprintf(' ... computing seeds')
        D = params.seedparams.weights;
        prc = params.seedparams.threshprc;
        thresh = percentile(D,prc);
        seeds = D > thresh;
        fprintf(' ... %g seeds selected',sum(seeds));
    end 
else
    fprintf(' ...using %g existing seeds', sum(params.seeds))
    seeds = params.seeds;
end

Ys = Y(:,seeds);
[m,n] = size(Ys);
fprintf(' => outer product with %g elements', n*n);
C = zeros(m,n*n);
% Loop through subjects
for i=1:m
    % Multiply seed region vector with its transpose
    Ci = Ys(i,:)'*Ys(i,:);
    %Ci = nk_PerfScale2(Ci,[],[],1);
    C(i,:) =Ci(:);
end

params.seeds = seeds;

end