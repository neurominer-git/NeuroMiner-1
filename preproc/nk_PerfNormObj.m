function [sY, IN] = nk_PerfNormObj(Y, IN)
% =========================================================================
% FORMAT function [sY, IN] = nk_PerfNormObj(Y, IN)
% =========================================================================
% Percentage scaling of single-subject data to the group mean. The group is
% indexed in IN.indY, where the row dimension indicates the group
% membership of the subjects and the col dimension the groups. If IN.indY
% is a vector with non-logical values a boolean group index meatrx is created.
% I/O arguments: 
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Nikolaos Koutsouleris, 10/2015

% =========================== WRAPPER FUNCTION ============================ 
if iscell(Y) && exist('IN','var') && ~isempty(IN)
    sY = cell(1,numel(Y)); 
    for i=1:numel(Y), 
        % Define active indices depending on training or testing situation
        if isfield(IN,'meanY')
           if isfield(IN,'TsInd'), IN.indY = IN.TsInd{i}; else IN.indY = []; end
        else
           if isfield(IN,'TrInd'), IN.indY = IN.TrInd{i}; else IN.indY = []; end
        end
        sY{i} =  PerfNormObj(Y{i}, IN ); 
    end
else
    % Define active indices depending on training or testing situation
    if isfield(IN,'meanY') 
       if isfield(IN,'TsInd'), IN.indY = IN.TsInd; else IN.indY = []; end
    else
       if isfield(IN,'TrInd'), IN.indY = IN.TrInd; else IN.indY = []; end
    end
    [ sY, IN ] = PerfNormObj(Y, IN );
end
% =========================================================================
function [sY, IN] = PerfNormObj(Y, IN)
global VERBOSE

% Defaults
% Zero-out non-finite features 
if ~isfield(IN,'zerooutflag') || isempty(IN.zerooutflag), IN.zerooutflag = 1;  end
if ~isfield(IN,'indY'), IN.indY = true(size(Y,1),1); end
if size(IN.indY,2) < 2 && ~islogical(IN.indY), 
    IN.indY = nk_MakeDummyVariables(IN.indY); 
end

nInd = find(any(IN.indY)); nG = numel(nInd); nX = size(IN.indY,2);

if ~isfield(IN,'meanY') || isempty(IN.meanY), 
    IN.meanY = zeros(nX,size(Y,2));
    indMeanYZero = false(nX,size(Y,2));
    meanexist=0; 
else
    meanexist=1; 
end

% Compute mean for current group as specified in indY(:,i)
% Check for values <= realmin in each group and remove these columns if
% all subjects of either group meet this criterion
if ~meanexist,
    for i=1:nG
        IN.meanY(nInd(i),:) = mean(Y(IN.indY(:,nInd(i)),:)) + 1;
        indMeanYZero(nInd(i),:) = IN.meanY(nInd(i),:) > realmin + 1;
    end
    IN.indL = sum(indMeanYZero(nInd,:)) == nG;
    if IN.indL
        if VERBOSE, fprintf(' Removed %g zero features in either of %g groups.',sum(~IN.indL), nG); end
    else
        IN.indL = true(1,size(Y,2));
    end
    
end

if ~meanexist, IN.meanY = IN.meanY(:,IN.indL); end
Y = Y(:,IN.indL);
sY = zeros(size(Y));

for i=1:nG
    if ~any(IN.meanY(nInd(i),:)) && nG == 1
        IN.meanY = mean(Y);
        sY = (Y*100./repmat(IN.meanY,size(Y,1),1))-100;
    else
        % Normalize single subject data of group to group mean
        sY(IN.indY(:,nInd(i)),:) = (Y(IN.indY(:,nInd(i)),:)*100./repmat(IN.meanY(nInd(i),:),sum(IN.indY(:,nInd(i))),1)) - 100;
    end
end

[ sY, IN ] = nk_PerfZeroOut(sY, IN);
