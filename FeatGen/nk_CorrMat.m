function [C, N, T, P, Pfdr ] = nk_CorrMat(Y,G,type)

[m1,n1] = size(Y);
[m2,n2] = size(G);
if m1~=m2
    error('Y and G must have same number of observations! Check your inputs!')
end
if ~exist('type','var')
    typ = 'pearson';
else
    typ = type;
end
T=[]; P=[];Pfdr=[];
switch typ
    case 'pearson'
        if ~isempty(which('fastcorr')); typ = 'pearson_fast'; end
    case 'spearman'
        %if ~isempty(which('spear')); typ = 'spearman_fast'; end
end
% Remove rows with Nan [if you don't want this, you may consider imputing
% data before running this script]

C = zeros(n2,n1); T = zeros(n2,n1); P = zeros(n2,n1); N = P;

for i=1:n2
    for j=1:n1
        
        tY = Y(:,j);
        tG = G(:,i);
        
        ix = ~isnan(tY) & ~isnan(tG);
        tY = tY(ix); tG = tG(ix);
        N(i,j) = sum(ix);
        
        switch typ
            case 'pearson'
                c = corrcoef(tY,tG);
                C(i,j) = c(1,2);
                if nargout > 1
                    [T(i,j), P(i,j)] = pvalPearson('b',C(i,j),m1);
                end                
            case 'pearson_fast'
                C(i,j) = fastcorr(tY,tG);
                if nargout > 1, [T(i,j), P(i,j)] = pvalPearson('b',C(i,j),m1); end
            case 'spearman'
                C(i,j) = corr(tY,tG,'type',typ);
            case 'spearman_fast'
                C(i,j) = spear(tY,tG);
        end   
    end
end

if nargout == 5
     [~,~,~, Pfdr] = fdr_bh(P,0.05);
end

function [t,p] = pvalPearson(tail, rho, n)
%PVALPEARSON Tail probability for Pearson's linear correlation.
t = rho.*sqrt((n-2)./(1-rho.^2)); % +/- Inf where rho == 1
switch tail
case 'b' % 'both or 'ne'
    try
        p = 2*spm_Tcdf(-abs(t),n-2);
    catch
        p= NaN;
    end
case 'r' % 'right' or 'gt'
    p = spm_Tcdf(-t,n-2);
case 'l' % 'left' or 'lt'
    p = spm_Tcdf(t,n-2);
end