function P = nk_TestPerfCompConnectivity( Y, L )

[m,n] = size(Y);
params.seedparams.threshprc = 95;
P=zeros(m,1);

% Perform LOO
for i = 1:m
   
    trainind = 1:m; trainind(i)=[];
    testind = i;
    fprintf('\nWorking on subject %g',i);
    
    %[Ytrain, meanY, stdY, meanY2] = nk_PerfStandardize2(Y(trainind,:),[],[],[],[],4);
    %test = nk_PerfStandardize2(Y(testind,:),[],meanY,stdY,meanY2,4);
    Ytrain = Y(trainind,:);
    Ytest = Y(testind,:);
    
    % Use FScore to weight voxels according to their predictive relevance
    D = nk_FScoreFeatRank(Ytrain,L(trainind));
    %mx = train_liblin(L, sparse(Y), '-c 2 -s 0');
    %D = mx.w;
    
    % Use weight vector to determine seed voxels and return outer product
    % of seed voxels (= structural covariance) of training population
    params.seedparams.weights = D;
    [ Ctrain, params ] = nk_PerfCompConnectivity( Ytrain, params );
    
    % Compute covariance matrix for held-back subject using seeds
    % determined in the training data (param.seeds)
    Ctest = nk_PerfCompConnectivity(Ytest,params);
    params = rmfield(params,'seeds');
    
    model = train_liblin(L(trainind),sparse(Ctrain),'-c 2 -s 2');
    
    P(i) = predict_liblin(L(testind),sparse(Ctest),model);
    
end  