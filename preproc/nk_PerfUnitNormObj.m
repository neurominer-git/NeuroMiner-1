function [sY, IN] = nk_PerfUnitNormObj(Y, IN)
% =========================================================================
% FORMAT function [sY, IN] = nk_PerfUnitNormObj(Y, IN)
% =========================================================================
% Normalized each feature in data matrix Y to a unit vector length
% I/O arguments:
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Nikolaos Koutsouleris, 10/2015

% =========================== WRAPPER FUNCTION ============================ 
if iscell(Y) && exist('IN','var') && ~isempty(IN)
   sY = cell(1,numel(Y)); 
    for i=1:numel(Y), [sY{i}, IN] = PerfUnitNormObj(Y{i}, IN); end
else
    [ sY, IN ] = PerfUnitNormObj(Y, IN );
end

% =========================================================================
function [sY, IN] = PerfUnitNormObj(Y, IN)

[mY,nY] = size(Y);
% Defaults
if isempty(IN),eIN=true; else eIN=false; end
% Zero-out non-finite features 
if eIN || ~isfield(IN,'zerooutflag')    || isempty(IN.zerooutflag), IN.zerooutflag = 2;  end

if eIN ||~isfield(IN,'normY') || isempty(IN.normY)
    if ~isfield(IN,'sIND') || isempty(IN.sIND),
        IN.sIND = true(mY,1); 
    else
        if ~islogical(IN.sIND), IN.sIND = logical(IN.sIND); end 
    end
    tY = Y(IN.sIND,:);
    if ~isfield(IN,'normMethod')  || isempty(IN.normMethod),  IN.normMethod = 1; end    
    if ~isfield(IN,'meanY')       || isempty(IN.meanY),       IN.meanY = mean(tY,1); end
    tY = bsxfun(@minus, tY, IN.meanY);
    IN.normY = zeros(1,nY);
    for i=1:nY
        IN.normY(i) = norm(tY(:,i),IN.normMethod);
    end
end
% Mean-center data
Y = bsxfun(@minus, Y, IN.meanY);
% Standardize matrix feature-wise to unit length
sY = bsxfun(@rdivide, Y,  IN.normY);
% Eliminate NaN features
[ sY, IN ] = nk_PerfZeroOut(sY, IN);