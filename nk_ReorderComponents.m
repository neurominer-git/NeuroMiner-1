function [cov_mat, reorder] = nk_ReorderComponents(T, S, method)
% T: the matrix to be reordered
% S: the source matrix which serves as template for the reordering
if ~exist('method','var') || isempty(method), method = 'mi'; end
switch method
    case 'mi'
        S = discretize(S);
        T = discretize(T);
end

switch method
    case 'cov'
        cov_mat = T'*S;
    otherwise
        si = size(S,2); ti = size(T,2);
        cov_mat = zeros(ti,si);
        for i=1:ti
            for j=1:si
                switch method
                    case 'mi'
                        cov_mat(i,j) = mi(T(:,i),S(:,j));
                    case 'pearson'
                        cov_mat(i,j) = fastcorr(T(:,i),S(:,j));
                    case 'spearman'
                        cov_mat(i,j) = corr(T(:,i),S(:,j),'type','spearman');
                end
            end
        end
end

[~,reorder] = max(cov_mat);
