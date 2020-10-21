function [sY, IN] = nk_PerfElimZeroObj(Y, IN)
% =========================================================================
% FORMAT function [sY, IN] = nk_PerfElimZeroObj(Y, IN)
% =========================================================================
% Remove features with zero-variance, and ANY Infs and NaNs.
% Furthermore, remove features with highly skewed distributions
% I/O arguments:
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Nikolaos Koutsouleris, 05/2018

% =========================== WRAPPER FUNCTION ============================ 
if ~exist('IN','var'), IN = []; end
if iscell(Y) 
    sY = cell(1,numel(Y)); 
    for i=1:numel(Y), [sY{i}, IN] =  PerfElimZeroObj(Y{i}, IN); end
else
    [ sY, IN ] = PerfElimZeroObj(Y, IN );
end

% =========================================================================
function [Y, IN] = PerfElimZeroObj(Y, IN)

global VERBOSE
% Defaults
if isempty(IN),eIN=true; else eIN=false; end

if eIN || ~isfield(IN,'NonPruneVec') || isempty(IN.NonPruneVec)
    if ~isfield(IN,'zero'), IN.zero = 1; end
    if ~isfield(IN,'inf'),  IN.inf = 1;  end
    if ~isfield(IN,'nan'),  IN.nan = 1;  end
    if ~isfield(IN,'perc'), IN.perc = []; end
    % Identify zero-variance columns
    if IN.zero == 1,  indNullVar = var(Y)==0; else, indNullVar = false(1,size(Y,2));end
    % Identify columns containing NaNs
    if IN.nan == 1,   indNan = any(isnan(Y)); else, indNan = false(1,size(Y,2)); end
    % Identify Inf columns
    if IN.inf == 1,   indInf = any(isinf(Y)); else, indInf = false(1,size(Y,2));end
    if ~isempty(IN.perc)
        R = nk_CountUniques(Y, IN.perc); indUT = R.UT';
    else
        indUT = false(1,size(Y,2)); 
    end
    % Put together
    IN.NonPruneVec = ~(indNullVar | indNan | indInf | indUT);
    if ~isempty(IN.perc)
        if VERBOSE,fprintf(' %g features eliminated (%g Zero-Var, %g NaNs, %g Infs, %g Single-value)', ...
        sum(~IN.NonPruneVec), sum(indNullVar), sum(indNan), sum(indInf), sum(indUT)); end
    else
        if VERBOSE,fprintf(' %g features eliminated (%g Zero-Var, %g NaN, %g Inf)', ...
        sum(~IN.NonPruneVec), sum(indNullVar), sum(indNan), sum(indInf)); end
    end
    
end
Y = Y(:,IN.NonPruneVec);