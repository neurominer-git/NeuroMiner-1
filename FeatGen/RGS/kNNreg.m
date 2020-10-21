function estY = kNNreg(X_test, X_train, Y_train, weights, extra_param, sorted_dists, NNs);

% Implementation of kNN regression algorithm, as used in
%
% A. Navot, L. Shpigelman, N. Tishby, E. Vaadia. Nearest Neighbor Based Feature Selection for Regression and its
% Application to Neural Activity. Submitted to NIPS 2005.
%
% estY = kNNreg(X_test, X_train, Y_train, weights, extra_param);
%
% input: X_train(i,j) is the value of feature j in training instance i.
%        Y_train(i) is the label of training instance i.
%        extra_param is a struct that may contain the following parameters for the algorithm:
%        k: number of neighbors (default: ceil(log2(# of instances)))
%        beta: beta is a Gaussian decay factor. (default as explained in the paper)
%
% output: estY(i) is the estimation of the function value for the i'th test instance.
%
%
% Can also be used in this format:
%
% estY = kNNreg(X_test, X_train, Y_train, weights, extra_param, sorted_dists, NNs);
%
% if the last two optional inputs are given, than the first two inputs are
% igonored. This mode is used by RGS to calculate the leave-one-out
% estimation. For this use, the these inputs should be:
%       sorted_dists(i,j) = the distance to from instance i in the training
%                           set the j'th closest instance in the training set 
%                          (excluding instance i itself)
%       NNs(i,j) = the index of the j'th closest instance to instance i in
%                  the training set (excluding instance i itself)
%       
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Written by Amir Navot & Lavi Shpigelman               
%% Date: June 3, 2005                        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isfield(extra_param, 'verbose'), extra_param.verbose = 0; end
if ~isfield(extra_param, 'k'), extra_param.k = ceil(log2(size(X_train,1))); disp(['k = ' num2str(extra_param.k)]); end
k = extra_param.k;
if ~isfield(extra_param, 'beta'),    
        %[sorted_dists, NNs, dists] = RGShelper(X_train, X_train.^2, ones(1,size(X_train,2)), 1:size(X_train,1));
        if nargin<7
            [sorted_dists, NNs] = getLooNN(X_train, X_train.^2, ones(1,size(X_train,2)), 1:size(X_train,1));
        end;
        extra_param.beta = mean(sorted_dists(:,floor(k/2)))/2;
        if extra_param.verbose, fprintf(' -> beta = %1.3f', extra_param.beta ); end
end

if nargin<7
  [sorted_dists, NNs] = findKNearestNeighborWeights(X_train, X_test, weights);
end
Y = Y_train(:);

estY = zeros(size(sorted_dists,1), 1);

coef = zeros(size(sorted_dists,1), k);
for ki = 1:k,
    coef(:,ki) = exp(-sorted_dists(:,ki) / extra_param.beta);
    estY = estY + Y(NNs(:,ki)) .* coef(:,ki);
end;
scoef = sum(coef,2);
estY = estY ./ scoef;



function [sorted_dists, NNs] = findKNearestNeighborWeights(X_train, X_test, weights);

P = X_train';
S = X_test';

alphas = weights;
alphas = alphas(:)';
Pnum = size(P,2);
Snum = size(S,2);

P2 =P.^2;
Pnorm = alphas * P2; % sum(P.^2,1);
PNorm = Pnorm(ones(1, Snum), :);
Snorm = (alphas * (S.^2))'; % sum(S.^2,1)';
SNorm = Snorm(:, ones(1,Pnum));

% Alphas = diag(alphas);
if Snum > 1,
    Alphas = spdiags(alphas',0,length(alphas),length(alphas));
    SP = 2 * (S'*Alphas) * P;
else,
    SP = 2 * (S'.*alphas) * P;
end


dists = SNorm + PNorm - SP;

[sorted_dists, NNs] = sort(dists,2);




