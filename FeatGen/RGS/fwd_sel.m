function selected_features = fwd_sel(X_train, Y_train, extra_param);


if ~isfield(extra_param, 'verbose'), extra_param.verbose=0;  end
if ~isfield(extra_param, 'k'), extra_param.k = ceil(log2(size(X_train,1))); if extra_param.verbose, disp(['k = ' num2str(extra_param.k)]); end; end
k = extra_param.k;
if ~isfield(extra_param, 'beta'),    
        [sorted_dists, NNs] = getLooNN(X_train, X_train.^2, ones(1,size(X_train,2)), 1:size(X_train,1));
        extra_param.beta = mean(sorted_dists(:,floor(k/2)))/2;
        if extra_param.verbose, disp(['beta = ' num2str(extra_param.beta )]); end
end
if ~isfield(extra_param, 'max_num_selected'), extra_param.max_num_selected = size(X_train,2); if extra_param.verbose, disp(['max_num_selected = ' num2str(extra_param.max_num_selected)]); end, end

num_feats = size(X_train,2);
m = size(X_train,1);
X_train_square = X_train.^2;

if (isfield(extra_param, 'initial_idx_feat'))
    idx_list = zeros(1,num_feats);
    idx_list(extra_param.initial_idx_feat) = 1;
    [sorted_dists, NNs, dists] = getLooNN(X_train, X_train_square, idx_list, 1:m);
    estY = kNNreg([], [], Y_train, idx_list, extra_param,sorted_dists, NNs);
    evalFunVal = -sum((Y_train-estY).^2);
    clear sorted_dists NNs;
else,
    idx_list = (zeros(1,num_feats)>0);
    dists = zeros(m);
    dists = dists + diag(realmax*ones(1,m));
end;

max_efv = -realmax;
updated=0;
for ns = 1:extra_param.max_num_selected,
    best_feat = 0;
    best_dists = zeros(size(dists));
    for feat = find(~idx_list),
        idx_list(feat)=~idx_list(feat);
        alt_dists = getAltDists(X_train, Y_train, dists, feat, idx_list(feat));
        [sorted_alt_dists, alt_NNs] = sort(alt_dists,2);
        estY = kNNreg([], [], Y_train, idx_list, extra_param,sorted_alt_dists, alt_NNs);
        alt_evalFunVal = -sum((Y_train-estY).^2);
        if alt_evalFunVal > max_efv,
            max_efv = alt_evalFunVal;
            best_feat = feat;
            best_dists = alt_dists;
            updated=1;
        end

        idx_list(feat)=~idx_list(feat);
        if extra_param.verbose && (mod(ns,1)==0)
        fprintf(['feat = ' num2str(feat) ' num_sel = ' num2str(ns) ' max_evf = ' ...
            num2str(max_efv) ' feats: ' num2str(find(idx_list)) '\n']);
        end

    end
    if updated==1
        idx_list(best_feat) = 1;
        dists = best_dists;
        updated=0;
    else
        if extra_param.verbose, disp('fwd_sel converged'); end
        break
    end

    if extra_param.verbose && (mod(ns,1)==0)
        fprintf(['num_sel = ' num2str(ns) ' max_evf = ' ...
            num2str(max_efv) ' feats: ' num2str(find(idx_list)) '\n']);
    end
end
evalFunVal = max_efv;
selected_features=find(idx_list);


function dists = getAltDists(X, Y, dists, changed_feat, added);

% dists = getAltDists(X, Y, dists, changed_feat, added);
% 
% input: X(i,j) is the value of feature j in instance i
%        Y(i) is the value of the ith example
%        dists is the current (square in two respects) distances
%        changed_feat is the index of the feature to flip
%        added is 1 if feature is added and 0 otherwise
%
%
% output:dists is the new distances matrix

P = X;
PL = Y;
S = P;
SL = PL;

m = size(X,1);
feat = X(:,changed_feat);
F1 = feat(:,ones(1,m));
F2 = F1';
delta = (F1 - F2).^2;
if added,
    dists = dists + delta;
else,
    dists = dists - delta;
end;
