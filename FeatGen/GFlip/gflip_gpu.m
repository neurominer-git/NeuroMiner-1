function [idx_feat, weights, score] = gflip_gpu(X_train, Y_train, extra_param)

% [idx_feat, score] = gflip(X_train, Y_train, extra_param);
%
% input: X_train(i,j) is the value of feature j in training instance i.
%        Y_train(i) is the label of training instance i.
%        extra_param is a struct that may contain the following parameters for the
%        algorithm:
%                  utility: 'linear', 'sigmoid' or 'zero-one'. (default is 'linear') 
%                  beta: if the utility is sigmoid, beta controles the utuility function u(m) = 1/(1+exp(-beta*m)). (default is beta=1)
%                  start_points: number of starting point to use in order to avoid local maxima. The 
%                                    runing time is linear in this number. (default is 1)
%                  max_iter: the number of pathes on all training data. (default is 1)
%                  initial_idx_list: idx_feat to start with. (default is empty list)
%                  verbose: 1 for verbose, 0 otherwise (default is 0)
%
% output: idx_feat contains the indecies of the features to be used
%         score is the value of the evaluation function for the features in
%         idx_feat

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Written by Amir Navot & Ran Gilad-Bachrach                  %%
%% Date: April 1, 2004                                          %%
%% Last update: November 4, 2004 by Ran Gilad-Bachrach          %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isfield(extra_param, 'verbose'), verbose = extra_param.verbose; else, verbose = 0;  end
if isfield(extra_param, 'start_points'), utility = extra_param.start_points; else, start_points = 1; if verbose,disp(['start_points = ' start_points]); end; end
best_score = realmin;
glabelsSet = GPUsingle(unique(Y_train));
gX_train = GPUsingle(X_train);
%gY_train = GPUsingle(Y_train);
for rep = 1:extra_param.start_points,
    
    rand('state',sum(100*clock));
    state = rand('state');
    
    [cur_weights, cur_idx_feat, cur_score] = gflipOneTime_gpu(gX_train, Y_train, extra_param, glabelsSet);
    if cur_score > best_score,
        best_score = cur_score;
        idx_feat = cur_idx_feat;
        weights = cur_weights;
    end;
end;
score = best_score;










