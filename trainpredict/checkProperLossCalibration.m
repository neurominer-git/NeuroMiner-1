function checkProperLossCalibration()

    N = 1e6;
    d = 3;
    x = (1/sqrt(d))*rand(N, d);

    %wTrue = rand(d,1);    
    %eta = 1./(1 + exp(-x * wTrue)); 

    eta = sqrt(sum(x - x.^2, 2));
    y = (rand(N,1) < eta);

    yhat = [];
    MODELS = {'linear' 'logistic'};

    for model = MODELS
        model = model{:};
        if strcmp(model, 'linear')
            w = ridge(y,x,0,0);        
            yhat = [yhat w(1) + x*w(2:end)];

        elseif strcmp(model, 'logistic')
            w = glmfit(x,y,'binomial');
            yhat = [yhat glmval(w,x,'logit')];
        end        
    end

    disp('comparison of mean y values:')
    disp([mean(y) mean(yhat)])    
    
    plotReliabilityDiagram(yhat, y, {'Linear' 'Logistic'}, 2);

    %figure;
    %plotCalibrationDeciles(yhat(:,2), y, 100);
