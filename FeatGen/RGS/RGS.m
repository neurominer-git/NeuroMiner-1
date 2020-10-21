function weights = RGS(X_train, Y_train, extra_param) 

% weights = RGS(X_train, Y_train, extra_param);
%
% Implementation of RGS algorithm, described in
%
% A. Navot, L. Shpigelman, N. Tishby, E. Vaadia. Nearest Neighbor Based Feature Selection for Regression and its
% Application to Neural Activity. Submitted to NIPS 2005.
%
% input: X_train(i,j) is the value of feature j in training instance i.
%        Y_train(i) is the label of training instance i.
%        extra_param is a struct that may contain the following parameters for the algorithm:
%             num_starts: number of starting point to use in order to avoid local maxima. The 
%                         runing time is linear in this number. (default is 2)
%
%             epochs: number of times to pass over all the training
%                     instances (default is 1)
%             k: number of neighbors (default: ceil(log2(# of instances)))
%             beta: beta is a Gaussian decay factor. (default as explained in the paper)
%             verbose: 1 for verbose, 0 otherwise (default is 0)
%
% output: weights(j) is the weight of the j's feature. higher is better.
%                  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Written by Amir Navot & Lavi Shpigelman 
%% Date: June 3, 2005
%% Last modification by Nikolaos Koutsouleris
%% Date: Oct 24, 2011
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global EVALFUNC

if ~isfield(extra_param, 'verbose'), extra_param.verbose = 0; end
if ~isfield(extra_param, 'num_starts'), extra_param.num_starts=2;  end;
if ~isfield(extra_param, 'k'), extra_param.k = ceil(log2(size(X_train,1))); end

%% Prepare
num_starts=extra_param.num_starts;
if isempty(EVALFUNC), EVALFUNC = 'MSE'; end
% Eliminate useless features
[ X_train, IN] = nk_PerfElimZeroObj(X_train, []);
fprintf('\n')
rand('state',sum(100*clock));
X_train_square = X_train.^2;
m=length(Y_train); 
switch EVALFUNC
    case {'MSE','RMSD','NRMSD','MAERR'}
        evalcrit = 'le';
    otherwise
        evalcrit = 'ge';
end
%% Print some info
fprintf('\nRGS algorithm paramaters:')
fprintf('\n+++++++++++++++++++++++++++++++++++++++++++++')
fprintf('\n# NN (k):\t\t%g', extra_param.k);
fprintf('\n# Starting points:\t%g', num_starts)
fprintf('\n# Epochs:\t\t%g', extra_param.epochs)
fprintf('\nTraining sample size:\t%g', numel(Y_train));
fprintf('\nHistogram analysis of target labels:')
[N,C] = hist(Y_train,10);
for i=1:numel(N)
    fprintf('\nBin center %g:\t%1.2f =>\t%g', i, C(i), N(i));
end
fprintf('\n+++++++++++++++++++++++++++++++++++++++++++++')
tic
for rep = 1:num_starts,
   flg = false;
   fprintf('\nIter %g/%g:', rep, num_starts)
   new_weights=RGSonce(X_train, Y_train, extra_param);
   [sorted_dists, NNs] = getLooNN(X_train, X_train_square, new_weights, 1:m);
   estY = kNNreg([], [], Y_train, new_weights, extra_param,sorted_dists, NNs);
   new_evalFunVal = feval(EVALFUNC, Y_train, estY); 
   if rep > 1, 
       tmpcrit = feval(evalcrit, new_evalFunVal, evalFunVal);
       if isnan(new_evalFunVal), continue, end
       if tmpcrit, flg = true; end
   else
       flg = true;
   end
   if flg
        evalFunVal = new_evalFunVal;
        weights = new_weights;
        fprintf(' ==> %s: %1.3f', EVALFUNC, new_evalFunVal)
   end
end
tElapsed = toc;
fprintf('\nDone in %1.2f seconds', tElapsed)
% Copy weights back to original unpruned space
tweights = zeros(numel(IN.NonPruneVec),1);
tweights(IN.NonPruneVec) = weights;
weights = tweights;
