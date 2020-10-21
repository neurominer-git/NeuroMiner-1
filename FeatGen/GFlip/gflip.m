function [idx_feat, weights, score] = gflip(X_train, Y_train, extra_param, cpus)

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

if ~isfield(extra_param, 'verbose'), extra_param.verbose = 0; end
if ~isfield(extra_param, 'start_points'), extra_param.start_points = 1; end
if strcmp(extra_param.utility,'sigmoid') && ...
    ((isfield(extra_param,'beta') && strcmp(extra_param.beta,'auto')) || ...
    ~isfield(extra_param,'beta')), 
    extra_param.beta = suggestBeta(X_train, Y_train);
    fprintf('\nAutodetected Beta=%g)',extra_param.beta)
elseif strcmp(extra_param.utility,'sigmoid')
    fprintf('\nPredefined Beta=%g)',extra_param.beta)
end

if nargin< 4
    try
        cpus = str2double(getenv('OMP_NUM_THREADS'));
    catch
        cpus =1;
    end
end
best_score = realmin;
labelsSet = unique(Y_train);
rep = 1;
while rep <= extra_param.start_points,
    
    rand('state',sum(100*clock));
    state = rand('state');
    
    [cur_weights, cur_idx_feat, cur_score] = gflipOneTime(X_train, Y_train, extra_param, cpus, labelsSet);
    if numel(cur_idx_feat) == 1,
        rep = rep - 1;
    end
    if cur_score > best_score,
        best_score = cur_score;
        idx_feat = cur_idx_feat;
        weights = cur_weights;
    end;
    rep = rep + 1;
end;
score = best_score;










