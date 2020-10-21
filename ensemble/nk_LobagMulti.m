function mED = nk_LobagMulti(E, T, L, C, uC, luC, G)

[dum, hE] = nk_MultiEnsPerf(E, T, L, C, G);
if nargin < 6
    uC = unique(C); luC = length(uC);
end

% Biased results
bias = hE ~= L;

% Biased variance
vb = 0; 

% Unbiased variance
vu = 0;

for i=1:luC
    vu = vu + var(E(~bias,C==uC(i)),1,2);     
    vb = vb + var(E(bias,C==uC(i)),1,2);
end
vu = vu / luC; vb = vb / luC;

if isempty(vb), vb = 0; end

% Unbiased variance
if isempty(vu), vu = 0; end
    
% Lobag computation
ED = mean(bias) + mean(vu) - mean(vb);

mED = mean(ED);
