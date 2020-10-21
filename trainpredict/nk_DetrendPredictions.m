function [Ycorr, beta, p, X] = nk_DetrendPredictions(beta, p, Y, X)

if exist('X','var') && ~isempty(X)
    E = Y - X; 
    % Compute mapping from prediction to its original label
    if ~exist('p','var') || isempty(p) 
        p = polyfit(X, Y, 1);
    end
    
    % Compute beta of rotation parameters between original label and error
    if ~exist('beta','var') || isempty(beta)
        IN.TrCovars = E;
        [~, IN] = nk_PartialCorrelationsObj(X, IN);
        beta = IN.beta;
    end
else
    IN.beta = beta;
end

% if original label does not exist, compute it using linear regression
% (polynomial coefficient p and predicted label are entered into equation)
if ~exist('X','var') 
    X = (Y - p(2)) / p(1); 
end
IN.TsCovars = X;
% Now detrend predictions based on X and beta
Ycorr = nk_PartialCorrelationsObj(Y, IN);

end