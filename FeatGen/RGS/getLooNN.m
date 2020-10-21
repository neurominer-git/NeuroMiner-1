function [sorted_dists, NNs, dists] = getLooNN(X, X_square, weights, inds)

% [sorted_dists, NNs] = getLooNN(X, X_square, weights, inds);
%
% is used by RGS and SKS
%
%       sorted_dists(i,j) = the distance to from instance i in the training
%                           set the j'th closest instance in the training set 
%                          (excluding instance i itself)
%       NNs(i,j) = the index of the j'th closest instance to instance i in
%                  the training set (excluding instance i itself)
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Written by Amir Navot & Lavi Shpigelman               
%% Date: June 3, 2005                        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

S = X(inds,:);
alphas = weights(:);
Pnum = size(X,1);
Snum = size(S,1);

Pnorm = X_square * alphas; % sum(P.^2,2);
PNorm = Pnorm(:,ones(1, Snum));
Snorm = ((S.^2) * alphas)'; % sum(S.^2,2)';
SNorm = Snorm(ones(1,Pnum), :);

if Snum > 1,
    %Alphas = diag(alphas);
    Alphas = spdiags(alphas,0,length(alphas),length(alphas));
    SP = 2 * X * ( Alphas * S' );
else
    SP = 2 * X * ( S' .* alphas );
end
dists = SNorm + PNorm - SP;
clear SNorm PNorm SP S ;
for k = 1:Snum, 
    dists(inds(k),k) = realmax;  % to avoid P_pos to be the same point
end;
[sorted_dists, NNs] = sort(dists',2);
if sum(sorted_dists<0)
  warning(['min(sorted_dists) = ' num2str(min(sorted_dists)) ' < 0']);
  sorted_dists = max(sorted_dists, 0);
%  error('sum(sorted_dists<0)')
end



