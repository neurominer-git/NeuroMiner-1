function [Ycorr, beta, p, X] = nk_DetrendPredictions2(beta, p, Y, X)
% =========================================================================
% function [Ycorr, beta, p, X] = nk_DetrendPredictions2(beta, p, Y, X)
% -------------------------------------------------------------------------
% INPUT
% beta : beta coefficients learned in a training sample
% p : polynomials
% X : Observed Label
% Y : Predicted Label
% OUTPUT
% Ycorr = adjusted Y
% =========================================================================
% (c) Nikolaos Koutsouleris, 01/2017
if exist('X','var') && ~isempty(X)
    
    E = Y - X; 
    % Compute mapping from prediction to its original label
    if ~exist('p','var') || isempty(p) 
        p = polyfit(X, Y, 1);
    end
    
    % Compute beta of rotation parameters between original label and error
    if ~exist('beta','var') || isempty(beta)
        IN.TrCovars = X;
        [Ecorr, IN] = nk_PartialCorrelationsObj(E, IN);
        beta = IN.beta;
    else
        IN.TrCovars = X;
        IN.TsCovars = X;
        IN.beta = beta;
        Ecorr = nk_PartialCorrelationsObj(E, IN);
    end
    
else
    % if original label does not exist, compute it using linear regression
    % (polynomial coefficient p and predicted label are entered into equation)
    if ~exist('X','var') 
        X = (Y - p(2)) / p(1); 
    end
    E = Y - X;
    IN.beta = beta;
    IN.TsCovars = X;
    Ecorr = nk_PartialCorrelationsObj(E, IN);
end

Ycorr = X + Ecorr; 

end