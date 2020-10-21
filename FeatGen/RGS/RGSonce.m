function weights = RGSonce(X_train, Y_train, extra_param)

m = size(X_train, 1); % Number of subjects
p = size(X_train, 2); % Number of dimensions

if ~isfield(extra_param, 'k'), extra_param.k = ceil(log2(size(X_train,1))); 
    if extra_param.verbose, fprintf(' k = %g', extra_param.k); end; end
k = extra_param.k;
if ~isfield(extra_param, 'epochs'), extra_param.epochs = 1; 
    if extra_param.verbose, fprintf(' epochs = %g', extra_param.epochs); end; end
if ~isfield(extra_param, 'beta'),    
        [sorted_dists, NNs] = getLooNN(X_train, X_train.^2, ones(1,size(X_train,2)), 1:size(X_train,1));
        extra_param.beta = mean(sorted_dists(:,floor(k/2)))/2;
        if extra_param.verbose, fprintf(' beta = %1.3f', extra_param.beta ); end
end

epochs = extra_param.epochs;
beta = extra_param.beta;

ws = ones(1, p);

X_train_square = X_train.^2;
t = 1;

for epoch = 1:epochs,
    
  perm = randperm(m);
  
  for pxi = 1:m,
    
    xi = perm(pxi);
 
    [sorted_dists, NNs] = getLooNN(X_train, X_train_square, ws.^2, xi);

    coef = zeros(k, 1);
    fhw = 0;
    du = zeros(1,p);
    dzw = zeros(1,p);
    for ki = 1:k,
      coef(ki) = exp(-sorted_dists(ki) / extra_param.beta);
      fhw = fhw + Y_train(NNs(ki)) .* coef(ki);
      du = du + ws .* (X_train(xi,:)- X_train(NNs(ki),:)).^2 * Y_train(NNs(ki)) * coef(ki);
      dzw = dzw + ws .* (X_train(xi,:)- X_train(NNs(ki),:)).^2 * coef(ki);
    end;
    scoef = sum(coef);
    u = fhw;
    fhw = fhw ./ scoef;
    zw = scoef;
    du = du *(-4/beta);
    dzw = dzw *(-4/beta);
    
    ws = ws +get_eta(t)*(Y_train(xi)-fhw)*(zw*du - u*dzw) / zw^2;
    t = t+1;
  end
end
weights = ws.^2;


function eta=get_eta(t);
eta = 1;