function yCI95 = nm_95confint(y)
N = size(y,1);                                      % Number of ‘Experiments’ In Data Set
%yMean = mean(y);                                   % Mean Of All Experiments At Each Value Of ‘x’
ySEM = nanstd(y)/sqrt(N);                              % Compute ‘Standard Error Of The Mean’ Of All Experiments At Each Value Of ‘x’
CI95 = tinv([0.025 0.975], N-1);                    % Calculate 95% Probability Intervals Of t-Distribution
yCI95 = bsxfun(@times, ySEM, CI95(:));              % Calculate 95% Confidence Intervals Of All Experiments At Each Value Of ‘x’