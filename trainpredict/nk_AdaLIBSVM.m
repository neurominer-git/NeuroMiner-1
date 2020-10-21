    N = length(X); % X training labels
    W = 1/N * ones(N,1); %Weights initialization
    M = 10; % Number of boosting iterations 

    for m=1:M

        %Calculate the error and alpha in adaBoost with cross validation
        cmd = ['-c ', num2str(C), ' -w ', num2str(W)];
        model = svmtrain(X, Y, cmd);
        [Xout, acc, ~] = svmpredict(X,Y,cmd);

        err = sum(.5 * W .* acc * N)/sum(W);
        alpha = log( (1-err)/err );

        % update the weight
        W = W.*exp( - alpha.*Xout.*X );
        W = W/norm(W);

    end