function [sY, IN] = nk_PerfExtDimObj(Y, IN)

% =========================== WRAPPER FUNCTION ============================ 
if iscell(Y) && exist('IN','var') && ~isempty(IN)
    sY = cell(1,numel(Y)); 
    for i=1:numel(Y), [sY{i}, IN] = PerfExtDimObj(Y{i}, IN); end
else
    [ sY, IN ] = PerfExtDimObj( Y, IN );
end

% =========================================================================
function [pY, IN] = PerfExtDimObj(Y, IN)

global VERBOSE

% Defaults
if isempty(IN),eIN=true; else eIN=false; end
[m,n] = size(Y);
% Check for and eliminate zero variance attributes
if eIN || ~isfield(IN,'indNonRem') || isempty(IN.indNonRem)
    if ~isfield(IN.mpp,'val') || IN.EXTDIM.PercMode == 1
        nD = 1:IN.opt;
    else
       switch IN.EXTDIM.PercMode
           case 2 % Percentage 
                nD = 1:floor( IN.opt * m );
           case 3 % Cumulative energy
                S = sum(IN.mpp.val)*IN.opt; 
                cS = cumsum(IN.mpp.val);
                nD = cS <= S ;
                if ~sum(nD), 
                    nD = false(1,numel(cS)); 
                    nD(1) = true;
                end
       end
    end
    IN.indNonRem = false(1,n); IN.indNonRem(nD) = true ; 
end
if VERBOSE, fprintf(' extracting %g components.', sum(IN.indNonRem)); end
try
    pY = Y(:,IN.indNonRem);
catch
    fprintf('problem')
end

