function [weights, score]  = nk_SimbaMain(X_train,  Y_train, extra_param, gpu)

% alphas = simba(X_train,  Y_train, extra_param);
%
% input: X_train(i,j) is the value of feature j in training instance i.
%               Y_train(i) is the label of training instance i.
%              extra_param is a struct that may contain the following parameters for the
%              algorithm:
%                  utility: 'linear' or 'sigmoid'. (default is 'linear') 
%                  beta: if the utility is sigmoid, beta controles the utuility function u(m) = 1/(1+exp(-beta*m)). (default is beta=1)
%                  start_points: number of starting point to use in order to avoid local maxima. The 
%                                    runing time is linear in this number. (default is 5)
%                  max_iter: the number of pathes on all training data. (default is 1)
%                  blocksize: the algorithm recaculate the distances after updating the features weights
%                             usig "blocksize" samples. default is 1. Use larger values in order to decrease the 
%                             running time.
%                  verbose: 1 for verbose, 0 otherwise (default is 0)
%
% output: weights(j) is the weight of the j's feature. higher is better.
%                  score is the value of the evaluation function for weights "weights"

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Written by Amir Navot & Ran Gilad-Bachrach                   %%
%% Date: April 1, 2004                                                                       %%
%% Last update: April 1, 2004 by Amir Navot                            %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isfield(extra_param, 'verbose'), verbose = extra_param.verbose; else, verbose = 0;  end
if isfield(extra_param, 'start_points'), start_points = extra_param.start_points; else, start_points = 5; if verbose, disp(['start_points = ' num2str(start_points)]); end ; end
if isfield(extra_param, 'utility'), utility = extra_param.utility; else, utility = 'linear'; if verbose, disp(['utility = ' utility]);end; end
if isfield(extra_param, 'beta'), 
    beta = extra_param.beta;
    if ischar(beta) && strcmp(beta,'auto'), beta = suggestBeta(X_train, Y_train); end
    extra_param.beta = beta;
else 
    beta = 1; if verbose disp(['beta = ' num2str(beta)]);end;
    extra_param.beta = beta;
end
if ~exist('gpu','var'), gpu = false; end

best_score = realmin;
weights = zeros(1,size(X_train,2));
X_train_square = X_train.^2;
m = size(X_train,1); ind = 1:m;
LabelSet = unique(Y_train);

for rep = 1:start_points,
    fprintf('\nRepetition %g/%g:',rep,start_points)
    rand('state',sum(100*clock));
    state = rand('state');
    if gpu        
        alphas = nk_SimbaLinearSigmoid_gpu(X_train,  Y_train, extra_param);
    else
        alphas = nk_SimbaLinearSigmoid(X_train,  Y_train, extra_param);
    end
    v = alphas.^2;
    v = v./max(v);
    [P_pos, P_neg, mindist_pos, mindist_neg, marginList] = nk_simbaHelper(X_train, X_train_square, Y_train, v, ind, LabelSet);
    
    if strcmp(utility,'linear')
        cur_score = abs(sum(marginList));
    elseif strcmp(utility,'sigmoid'),
        cur_score = sum(1./(1 + exp(-beta * marginList)));
    else
        error('Illegal utility function, should be linear or sigmoid');
    end
    if cur_score > best_score,
        best_score = cur_score;
        weights = v';
    end
    fprintf('\nBest score after %g repetition(s): %1.2f',rep, best_score) 
end
score = best_score;
