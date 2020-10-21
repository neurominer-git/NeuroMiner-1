function CI = polyparci(PolyPrms,PolyS,alpha)
% POLYPARCI takes PolyPrms, the 1xN vector of parameter estimates from polyfit,
%   PolyS, the structure returned by polyfit, and alpha, and returns the
%   alpha confidence intervals for the parameters in CI.  The default
%   for alpha is 0.95.  
% 
%   RETURNS: CI (2xN matrix of confidence intervals whose columns match the 
%            1xN vector of coefficients returned by polyfit)
% 
% CALLING POLYPARCI: 
% [p,S,mu] = polyfit(x,y,n);    % Fit the data (polyparci does not use mu)
% ci = polyparci(p,S);          % NOTE: specifying a value for alpha is 
%                                 optional
% 
% Star Strider  2012 07 05;    
% Update         2014 02 10 (Compatible with 2013b changes. Corrected typo. 
%                   Inverted CI matrix to 2xN be compatible with 1xN row 
%                   vector of coefficients polyfit produces.)
% 

% Check for a specified value of alpha: 
if nargin < 3
    alpha = 0.95;
end

% Check for out-of-range values for alpha and substitute if necessary: 
if alpha < 1.0E-010
    alpha = 1.0E-010;
elseif alpha > (1 - 1.0E-010)
    alpha = 1 - 1.0E-010;
end

% Calculate the covariance matrix of the parameters (COVB) and the standard
%   errors (SE) from the S structure.  (See the polyfit documentation.)   

COVB = (PolyS.R'*PolyS.R)\eye(size(PolyS.R)) * PolyS.normr^2/PolyS.df;

SE = sqrt(diag(COVB));                              % Standard Errors

[PrmSizR,PrmSizC] = size(PolyPrms);                 % Convert parameter vector to column if necessary
if PrmSizR < PrmSizC
    PolyPrms = PolyPrms';
end

% Calculate t-statistic and inverse t-statistic:
%       tstat is an implicit function that calculates the probability
%       from the cumulative t-distribution t_cdf for various values of
%       tval supplied to it by the fzero call until it converges to the
%       required value of alpha.  The t_cdf function calculates the
%       cumulative t-distribution.  The two lines that calculate tstat
%       and T together calculate the inverse t-distribution, returning it
%       as T.

tstat = @(tval) (alpha - t_cdf(tval,PolyS.df) );    % Function to calculate t-statistic for p = alpha and v = PolyS.df

[T,fval] = fzero(tstat, 1);                         % Calculate t-statistic for p = alpha and v = PolyS.df

T = abs(T);                                         % Critical +ve value from t-distribution

ts = T * [-1  1];                                   % Create vector of  t-statistic values

CI  = bsxfun(@plus,bsxfun(@times,SE,ts),PolyPrms)'; % Confidence Intervals Matrix (2xN)


% CALCULATE THE CUMULATIVE T-DISTRIBUTION: 
    function PT = t_cdf(t,v)
        % t_cdf(t,v) calculates the cumulative t-distribution probability
        %   given the t-statistic t and degrees-of-freedom v.  The
        %   routine to calculate the inverse t-distribution uses this,
        %   tstat, and the fzero call.  Compared to the Statistics
        %   Toolbox function tcdf and tinv, t_cdf and T have
        %   relative errors of about 1E-12.
        
        IBx = v./(t.^2 + v);                      % x for IxZW  NOTE: function of t-statistic
        IBZ = v/2;                                % Z for IxZW
        IBW = 0.5;                                % W for IxZW
        Ixzw = betainc(IBx, IBZ, IBW);            % Incomplete beta function
        PT = 1-0.5*Ixzw;                          % Cumulative t-distribution        
    end


end

% =========================  END: polyparci.m  =========================