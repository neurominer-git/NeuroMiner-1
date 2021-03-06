function [weights, idx_feat, score] = gflipOneTime_gpu(gX_train, Y_train, extra_param, glabelsSet)

% [idx_feat, score] = gflip_incdists(obj1, X_train, Y_train, extra_param);
%
% input: X_train(i,j) is the value of feature j in training instance i.
%        Y_train(i) is the label of training instance i.
%        extra_param is a struct that may contain the following parameters for the
%        algorithm:
%                  utility: 'linear', 'sigmoid' or 'zero-one'. (default is 'linear') 
%                  beta: if the utility is sigmoid, beta controles the utuility function u(m) = 1/(1+exp(-beta*m)). (default is beta=1)
%                  initial_idx_list: idx_feat to start with. (default is empty list)
%                  max_iter: the number of pathes on all training data. (default is 1)
%                  verbose: 1 for verbose, 0 otherwise (default is 0)
%
% output: idx_feat contains the indecies of the features to be used
%         score is the value of the evaluation function for the features in
%         idx_feat

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Written by Amir Navot & Ran Gilad-Bachrach                   %%
%% created: May 5, 2004                                          %%
%% Last update: November 4, 2004 by Amir Navot                   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isfield(extra_param, 'verbose'), verbose = extra_param.verbose; else, verbose = 0; end
if isfield(extra_param, 'utility'), utility = extra_param.utility; else, utility = 'linear'; if verbose, disp(['utility = ' utility]); end; end
if isfield(extra_param, 'beta'), beta = extra_param.beta; else, beta = 1; if verbose, disp(['beta = ' num2str(beta)]); end; end
if isfield(extra_param, 'max_iter'), max_iter = extra_param.max_iter; else, max_iter = 1; if verbose, disp(['max_iter = ' num2str(max_iter)]); end; end
if isfield(extra_param, 'beta'), beta = extra_param.beta; else, beta = 1; if verbose, disp(['beta = ' num2str(beta)]); end; end

feat_num = size(gX_train,2);
weights = zeros(1,feat_num);
m = size(gX_train,1);
gX_train_square = gX_train.^2;
%  idx_list = (rand(1,feat_num)>0.5);

if (isfield(extra_param, 'initial_idx_feat'))
    idx_list = zeros(1,feat_num);
    idx_list(extra_param.initial_idx_feat) = 1;
    [P_pos, P_neg, mindist_pos, mindist_neg, marginList,dists] = simbaHelperGpu(gX_train, gX_train_square, Y_train,  idx_list, 1:m);
    if strcmp(utility,'linear'),
        score = sum(marginList);
    elseif strcmp(utility,'sigmoid'),
        score = sum(1./(1 + exp(-beta * marginList)));
    elseif strcmp(utility,'zero-one'),
        score = sum(0.5*(sign(marginList)+1));
    else
        error('Illegal utility function, should be linear or sigmoid');
    end
else
    idx_list = (zeros(1,feat_num)>0);
    score = 0;
    dists = zeros(m);
    dists = dists + diag(realmax*ones(1,m));
end;
%suggestedBeta = zeros(max_iter,1);
%figure(10);clf;hold on
dists = GPUsingle(dists);
for iter = 1:max_iter,
    tic;
    perm = randperm(feat_num);
    updated = 0;
    
    for ii = 1:feat_num,
        
        feat = perm(ii);
        idx_list(feat)=~idx_list(feat);
        
        [P_pos, P_neg, mindist_pos, mindist_neg, marginList, alt_dists] = ...
            gflipHelperGpu(gX_train, Y_train, dists, feat, idx_list(feat), glabelsSet);
        
        if strcmp(utility,'linear'),
            alt_score = sum(marginList);
        elseif strcmp(utility,'sigmoid'),
            alt_score = sum(1./(1 + exp(-beta * marginList)));
        elseif strcmp(utility,'zero-one'),
            alt_score = sum(0.5*(sign(marginList)+1));
        else
            error('Illegal utility function, should be linear or sigmoid');
        end
        
        weights(feat) = (alt_score - score) * (2*idx_list(feat)-1);
        if alt_score > score,
            
            score = alt_score;
            updated=1;
            dists = alt_dists;
        else
            idx_list(feat)=~idx_list(feat);
        end;
       
        %      plot(1:length(idx_list), idx_list,'*');
        %if verbose && (mod(ii,10)==0)
        %    fprintf(['epoch = ' num2str(iter) ', iter= ' num2str(ii) ', score= ' ...
        %            num2str(score) ' active ' num2str(sum(idx_list)) '\n']); 
        %end
    end
    tElapsed=toc;
    if updated==0,
        fprintf(' converged in %g iterations. ',iter);
        break;
    else
        fprintf('\nScore: %g ... in %1.2f sec', score, tElapsed);
    end;
    %plot(ii,sort(weights),'b*'); drawnow;
end
idx_feat=find(idx_list);
%mSuggestedBeta = mean(suggestedBeta);
if ~( any(idx_list>0 ~= weights>0) & any(idx_list>0 ~= weights>=0) ),
    warning('G-flip - Internal error');
end





