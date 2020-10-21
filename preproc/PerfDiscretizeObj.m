function [ dY, IN ] = PerfDiscretizeObj(Y, IN)

% Defaults
if isempty(IN),eIN=true; else eIN=false; end
% Default params for binning
if eIN|| ~isfield(IN,'DISCRET') || isempty(IN.DISCRET),  
    IN.DISCRET.binstart = 0;
    IN.DISCRET.binsteps  = 0.5;
    IN.DISCRET.binstop  = 5;
end
% Compute mean if not in IN.mY
if eIN || ~isfield(IN,'mY') || isempty(IN.mY),
    IN.mY = mean(Y); 
end
% Compute mean if not in IN.sY
if eIN || ~isfield(IN,'sY') || isempty(IN.sY),  
    IN.sY = std(Y);
end

[m,n] = size(Y);
dY =zeros(m,n);

alphas = IN.DISCRET.binstart : IN.DISCRET.binsteps : IN.DISCRET.binstop;

tmY = repmat(IN.mY,m,1); tsY = repmat(IN.sY,m,1);

for j=1:numel(alphas)-1
    dY( Y > ( tmY + alphas(j) .* tsY) ) = alphas(j+1);
    dY( Y < ( tmY - alphas(j) .* tsY) ) = -alphas(j+1);
end

dY( Y > tmY + alphas(end) .* tsY ) = IN.DISCRET.binstop;
dY( Y < tmY - alphas(end) .* tsY ) = -(IN.DISCRET.binstop);

