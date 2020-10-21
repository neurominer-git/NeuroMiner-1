function IN = nk_PerfANOVAObj(Y, IN)
global VERBOSE
[~,n] = size(Y);
IN.p = zeros(n,1);
IN.F = zeros(n,1);
IN.R2 = zeros(n,1);
if VERBOSE, fprintf('\t running ANOVA on %g variables ',n); end

for i=1:size(Y,2)
    
    if VERBOSE;fprintf('.'); end
    try
        IN.beta = pinv(IN.X)*Y(:,i); % almost as above
    catch
        fprintf('error')
    end
    % here we use a pseudoinverse:
    % X is rank deficient, i.e. regressors are not independent, since any linear combination of
    % 3 columns can give us the 4th one, thus X'X is also rank deficient and singular ie inv(X'*X)
    % doesn't exist -- there is no matrix A such as A(X'X) = I - however pinv gives a unique
    % solution that minimizes the square distance to the data.
    % see
    % http://mathworld.wolfram.com/Moore-PenroseMatrixInverse.html
    % http://en.wikipedia.org/wiki/Moore_Penrose_pseudoinverse

    IN.Yhat = IN.X*IN.beta;
    IN.df = rank(IN.X)-1;
    

    IN.Res = Y(:,i) - IN.Yhat;
    IN.SStotal = norm(Y(:,i) - mean(Y(:,i))).^2;
    IN.SSeffect = norm(IN.Yhat - mean(IN.Yhat)).^2;
    IN.SSerror  = norm(IN.Res - mean(IN.Res)).^2;
    IN.dferror = length(Y(:,i)) - IN.df - 1;
    IN.R2(i) = IN.SSeffect / IN.SStotal;
    IN.F(i) = (IN.SSeffect / IN.df) / (IN.SSerror / IN.dferror);
    IN.p(i) = 1 - fcdf(IN.F(i),IN.df,IN.dferror);

end
