function [sY, IN] = nk_PerfRemMeanDiffObjOld(Y, IN)
% =========================================================================
% FORMAT function [sY, IN] = nk_PerfRemMeanDiffObj(Y, IN)
% =========================================================================
% Normalizes given data (specified ind dIND) to global mean of given groups 
% (specified in sIND) by first subtracting global mean (meanY) and then 
% subtracting offsets (meanG) from data 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Nikolaos Koutsouleris, 10/2015

% =========================== WRAPPER FUNCTION ============================ 
if iscell(Y) && exist('IN','var') && ~isempty(IN)
    sY = cell(1,numel(Y)); 
    for i=1:numel(Y), 
        % Define active indices depending on training or testing situation
        if isfield(IN,'meanY') && isfield(IN,'meanG')
           if isfield(IN,'sTsInd'), IN.sIND = IN.sTsInd{i}; else IN.sIND =[]; end
           if isfield(IN,'dTsInd'), IN.dIND = IN.dTsInd{i}; else IN.dIND =[]; end
        else
           if isfield(IN,'sTrInd'), IN.sIND = IN.sTrInd{i}; else IN.sIND = []; end
           if isfield(IN,'dTrInd'), IN.dIND = IN.dTrInd{i}; else IN.dIND = []; end
        end
        [ sY{i}, IN ] = PerfRemMeanDiffObj(Y{i}, IN ); 
    end
else
    % Define active indices depending on training or testing situation
    if isfield(IN,'meanY') && isfield(IN,'meanG')
        if isfield(IN,'sTsInd'), IN.sIND = IN.sTsInd; else IN.sIND =[]; end
        if isfield(IN,'dTsInd'), IN.dIND = IN.dTsInd; else IN.dIND =[]; end
    else
        if isfield(IN,'sTrInd'), IN.sIND = IN.sTrInd; else IN.sIND = []; end
        if isfield(IN,'dTrInd'), IN.dIND = IN.dTrInd; else IN.dIND = []; end
    end
    [ sY, IN ] = PerfRemMeanDiffObj(Y, IN );
end

% =========================================================================
function [sY, IN] = PerfRemMeanDiffObj(Y, IN)
% Defaults
global VERBOSE
if isempty(IN),eIN=true; else eIN=false; end

if eIN || ~isfield(IN,'sIND') || isempty(IN.sIND), IN.sIND = true(size(Y,1),1); end
if ~islogical(IN.sIND) 
    if size(IN.sIND,2)<2, 
        IN.sIND = nk_MakeDummyVariables(IN.sIND); 
    else
        IN.sIND = logical(IN.sIND);
    end
end
% Index vector / logicals matrix of subjects to apply the correction
% parameters to.
if eIN || ~isfield(IN,'dIND') || isempty(IN.dIND), IN.dIND = true(size(Y,1),1); end
if ~islogical(IN.dIND) 
    if size(IN.dIND,2)<2, 
        IN.dIND = nk_MakeDummyVariables(IN.dIND); 
    else
        IN.dIND = logical(IN.dIND);
    end
end

nGM = size(IN.dIND,2); if nGM == 1, IN.dIND = [IN.dIND ~IN.dIND]; nGM = 2; end
nGY = size(IN.sIND,2); if nGY == 1, IN.sIND = [IN.sIND ~IN.sIND]; nGY = 2; end

if ~isfield(IN,'meanY') || ~isfield(IN,'meanG') || isempty(IN.meanY) || isempty(IN.meanG)
    
    indG = any(IN.sIND,2);
    [~, D] = size(Y(indG,:));
    
    % Compute overall mean of the data
    IN.meanY = mean(Y(indG,:));

    %sYx = Y(indG,:) - repmat(meanY,S,1);
    
    % Compute group-specific means
    IN.meanG = zeros(nGM,D);
    for i = 1:nGM
        if ~any(IN.dIND(:,i)), continue, end
        IN.meanG(i,:) = mean(Y(IN.sIND(indG,i),:)); 
    end

end

% Mean centering of data to global mean
%fprintf('\nRemoving global mean')
%sY = Y - repmat(meanY, Sa, 1);
sY=Y;
onceflag = false;

for i=1:nGM
    dINDi = IN.dIND(:,i);
    dINDi2 = IN.dIND; dINDi2(:,i)=[]; LGO = ~any(dINDi2(:));
    if ~any(IN.dIND(:)) 
        if VERBOSE, fprintf('\nNo data mean offset found for these data. Computing it right now!'); end
        % This is data whose group mean has not been computed from training data, 
        % because training data was not available ==> OOCV data
        % in this case compute offset right now
        meanGi = mean(sY) - IN.meanY;
        dINDi = true(size(sY,1),1);
        onceflag = true;
    elseif ~any(dINDi) 
        if VERBOSE, fprintf('\nSkipping offset correction because of no observation in group %g.',i); end
        continue,
    elseif LGO
        if VERBOSE, fprintf('\nApplying offset correction to leave-group-out data %g.',i); end
        meanGi = mean(sY) - IN.meanY;
    else
        meanGi = IN.meanG(i,:) - IN.meanY;
    end
    
    % Normalize entire data to global mean of group i
    sY(dINDi,:) = bsxfun(@minus, sY(dINDi,:), meanGi);
    
    if onceflag, break, end;
    
end

sY(~isfinite(sY))=0;