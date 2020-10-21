% Bins probability estimates and also returns the empirical probability
% estimate within each bin
function [pDisc, empProb] = binProbabilities(p, y, gridSize)
    
    % Sort probability estimates in increasing order
    p = reshape(p, [1 numel(p)]);
    [p, i] = sort(p, 'ascend');
    y = y(i);
    
    pDisc = round(10^gridSize * p)/10^gridSize; % Discretize p to 2dp, to perform binning    
    
    k = max(y); % k is the value of the label we are checking calibration against
    empProb = zeros(size(pDisc)); % The values plotted on the y-axis, empirical probability estimates

    for p0 = unique(pDisc)
        % For all the p values equal to p0, what's the empirical frequency
        % of the event yi = k?
        idx = pDisc == p0;
        empProb(idx) = mean(y(idx) == k);
    end

    % Revert to original indices (no longer in sorted probability order)
    [~,j] = sort(i);
    pDisc = pDisc(j);
    empProb = empProb(j);

    pDisc = reshape(pDisc, [numel(pDisc) 1]);
    empProb = reshape(empProb, [numel(empProb) 1]);
