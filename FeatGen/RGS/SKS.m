function scores = SKS(X_train, Y_train, extra_param);

% scoress = SKS(X_train, Y_train, extra_param);
%
% Implementation of SKS algorithm, described in
%
% A. Navot, L. Shpigelman, N. Tishby, E. Vaadia. Nearest Neighbor Based Feature Selection for Regression and its
% Application to Neural Activity. Submitted to NIPS 2005.
%
% input: X_train(i,j) is the value of feature j in training instance i.
%        Y_train(i) is the label of training instance i.
%        extra_param is a struct that may contain the following parameters for the algorithm:
%             k: number of neighbors (default: ceil(log2(# of instances)))
%             beta: beta is a Gaussian decay factor. (default as explained in the paper)
%
% output: scores(j) is the weight of the j's feature. higher is better.
%                  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Written by Amir Navot & Lavi Shpigelman               
%% Date: June 3, 2005                        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isfield(extra_param, 'verbose'), extra_param.verbose=0;  end
if ~isfield(extra_param, 'k'), extra_param.k = ceil(log2(size(X_train,1))); if extra_param.verbose, disp(['k = ' num2str(extra_param.k)]); end; end
k = extra_param.k;
if ~isfield(extra_param, 'beta'),    
        [sorted_dists, NNs] = getLooNN(X_train, X_train.^2, ones(1,size(X_train,2)), 1:size(X_train,1));
        extra_param.beta = mean(sorted_dists(:,floor(k/2)))/2;
        if extra_param.verbose, disp(['beta = ' num2str(extra_param.beta )]); end
end

m = size(X_train, 1);
p = size(X_train, 2);
scores = zeros(1, p);


[s I]=sort(X_train,1);
curNN_dists=zeros(k);
curNNs=zeros(k);

% intialization
for pi=1:p
  last_nn(pi)=k+1;
  curNNs=I(2:k+1,pi);
  curNN_dists=(s(2:k+1,pi)-s(1,pi)).^2;
  first_nn(pi)=1;
  ys = Y_train(I(1,pi));
  
  pred = sum(Y_train(curNNs) .* exp(-curNN_dists/extra_param.beta)) ./ sum(exp(-curNN_dists/extra_param.beta));
  curw(pi) = (pred - ys).^2;

end
scores = scores - curw;
for mi=2:m
  for pi=1:p
    if last_nn(pi)==m
      if first_nn(pi)>mi | last_nn(pi)<mi
        error('first_nn(pi)>mi | last_nn(pi)<mi')

      elseif first_nn(pi) == mi
        range=first_nn(pi)+1:last_nn(pi);
      else
        range=[first_nn(pi):mi-1 mi+1:last_nn(pi)];
      end

    else
      while last_nn(pi)<m & s(last_nn(pi)+1,pi)-s(mi,pi)<s(mi,pi)-s(first_nn(pi),pi)
        first_nn(pi)=first_nn(pi)+1;
        last_nn(pi)=last_nn(pi)+1;
      end
    end
    if last_nn(pi)<mi
      first_nn(pi)=first_nn(pi)+1;
      last_nn(pi)=last_nn(pi)+1;
      range=first_nn(pi):last_nn(pi)-1;
    elseif first_nn(pi) == mi
      range=first_nn(pi)+1:last_nn(pi);
    else
      range=[first_nn(pi):mi-1 mi+1:last_nn(pi)];
    end
    curNNs=I(range,pi);
    curNN_dists=(s(range,pi)-s(mi,pi)).^2;
    ys = Y_train(I(mi,pi));
    pred = sum(Y_train(curNNs) .* exp(-curNN_dists/extra_param.beta)) ./ sum(exp(-curNN_dists/extra_param.beta));
    curw(pi) = (pred - ys).^2;

  end
  scores = scores - curw;
end


