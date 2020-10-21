function [sY, mI, mIG, H] = nk_MultiCentIntensNorm(Y, V, G, B, mI, mIG, flags)
% This function performs variable normalization for multi-center studies
% Inputs:
% Y = data matrix for offset correction
% V = overlap search variable
% G = Group index vector
% B = binwidth for histogram analysis of V

if ~exist('flags','var') || isempty(flags), 
    flags.display = true; 
    flags.meancenter = true;
end
if ~exist('B','var') || isempty(B), B = 10; end
if ~exist('G','var') || isempty(G), error('Group vector missing or undefined'); end
if ~exist('V','var') || isempty(V), error('Overlap search vector missing or undefined'); end
if ~exist('Y','var') || isempty(Y), error('Data matrix missing or undefined'); end
if ( ~exist('mI','var') || isempty(mI) ) || ~exist('mIG','var') || isempty(mIG),
    flags.train = true;
else
    flags.train = false; 
end

if ndims(G) == 1, % Check if G is a single vector with number indiciating centers / groups
   U = unique(G); lx = length(G); nclass = length(U); 
   tO = false(lx,nclass); tO(sub2ind([lx nclass], 1:lx, G')) = true;
   G = tO;
elseif ~islogical(G)
    G = logical(G);
end

nG = size(G,2);

if flags.train
    
    %% Overlay variable histograms and find bin with maximum sample size 
    N = zeros(B,nG); H = zeros(B,nG); C = zeros(B,nG);
    Vvec = min(V) : ceil((max(V) - min(V))/B) : max(V);

    for i = 1 : nG
       [N(:,i),  C(:,i)] = hist(V(G(:,i)),B) ;
       H(:,i) = N(:,i)*100 / sum(G(:,i)); 
    end
    
    %% Compute cumulative histogram across groups / centers
    sH = sum(N,2);
    [sH_max, sH_max_ind] = max(sH);
    fprintf('\nMaximum number of cases detected in bin #%g (%g).',sH_max_ind,sH_max)

    if flags.display, figure;bar(sH); end

    %% Compute group offsets in this bin
    if sH_max_ind == B
        I = V <= Vvec(end) & V > Vvec(end-1);
    else
        I = V >= Vvec(sH_max_ind) & V < Vvec(sH_max_ind+1);
    end
    
    %% Compute overall mean in this bin
    mI = mean(Y(I,:));
    mIG = zeros(nG, size(Y,2));
    
    %% Eventually, mean center entire data
    if flags.meancenter, 
        Y = bsxfun(@minus, Y, mI); 
        %% Compute group / center offsets from globally mean centered data
        for i = 1 : nG
            ind = I & G(:,i);
            mIG(i,:) = mean(Y(ind,:));
        end
    else
        %% Compute group / center offsets from data mean
        for i = 1 : nG
            ind = I & G(:,i);
            mIG(i,:) = mean(Y(ind,:)) - mI;
        end
    end
else
    Y = bsxfun(@minus, Y, mI); 
end

%% Now subtract the offsets from the entire groups
sY = zeros(size(Y));
for i = 1 : nG
    ind = G(:,i); 
    sY(ind,:) = Y(ind,:) - repmat(mIG(i,:),sum(ind),1); 
end

end