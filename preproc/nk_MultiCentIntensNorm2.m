function [sY, params] = nk_MultiCentIntensNorm2(Y, V, G, B, Brange, IN, flags)
% This function performs variable normalization for multi-center studies
% Inputs:
% Y = data matrix for offset correction
% V = overlap search variable
% G = Group index vector
% B = binwidth for histogram analysis of V

if ~exist('flags','var') || isempty(flags), 
    flags.train = 1;
    flags.display = true; 
    flags.corrmethod = 1;
end
if ~exist('B','var') || isempty(B), B = 10; end
if ~exist('G','var') || isempty(G), error('Group vector missing or undefined'); end
if ~exist('V','var') || isempty(V), error('Overlap search vector missing or undefined'); end
if ~exist('Y','var') || isempty(Y), error('Data matrix missing or undefined'); end

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
    
    % Build edges vector for histc
    fprintf('\n\n*********************************************************************')
    fprintf('\n*** Correcting group offset by finding maximum group overlap in V ***')
    fprintf('\n*********************************************************************')
    if exist('Brange','var') && ~isempty(Brange) && Brange(1) >= min(V) && Brange(2) <= max(V)
        Edges = Brange(1) : (Brange(2) - Brange(1))/B : Brange(2);
    else
        Edges = min(V) : (max(V) - min(V))/B : max(V);
    end
    nE = numel(Edges);
    N = zeros(nE,nG); 
    H = zeros(nE,nG); 
    % Compute abolute num (N) and relative num (H) histograms for each center
    for i = 1 : nG
       N(:,i) = histc( V( G(:,i) ), Edges) ;
       H(:,i) = N(:,i) * 100 / sum( G(:,i) ); 
    end
    fprintf('\n\nEdges for histogram computation')
    fprintf('\n%g', Edges)
    
    %% Compute cumulative histogram across groups / centers
    fprintf('\n\nRelative number histogram:')
    for i =1 : nE
        fprintf('\nBin %g\t=>',Edges(i))
        fprintf('%1.0f\t',H(i,:))
    end
    fprintf('\n\nAbsolute number histogram:')
    for i =1 : nE
        fprintf('\nBin %g\t=>',Edges(i))
        fprintf('%1.0f\t',N(i,:))
    end
    sH = sum(H,2); [sH_max, sH_max_ind] = max(sH);
    fprintf('\n\nMaximum number of cases (N = %g) detected in bin #%g (%g).', sum(N(sH_max_ind,:),2), sH_max_ind,sH_max)
    %% Compute index to subjects used to train offset model
    if sH_max_ind == nE,
        I = V >= Edges(sH_max_ind-1);
    else
        I = V >= Edges(sH_max_ind) & V < Edges(sH_max_ind+1);
    end
    
    switch flags.corrmethod 
        case 1
            fprintf('\nCorrecting centers offsets using Partial correlations')
            IN.TrCovars = G(I,:); [~, IN] = nk_PartialCorrelationsObj(Y(I,:),IN);
            params = struct('Edges',Edges,'H',H,'N',N,'SelEdgeIndex',sH_max_ind,'I',I,'beta',IN.beta);
        case 2
            fprintf('\nCorrecting centers offsets using normalization of center data to each center''s mean value')
            [~, meanY, indL] = nk_PerfNorm(Y(I,:), G(I,:));
            params = struct('Edges',Edges,'H',H,'N',N,'SelEdgeIndex',sH_max_ind,'I',I,'meanY',meanY,'indL',indL);
        case 3
            fprintf('\nCorrecting centers offsets using center offset correction to global mean')
            [~, meanY, meanG] = nk_PerfRemMeanDiff(Y(I,:), G(I,:), G(I,:));
            params = struct('Edges',Edges,'H',H,'N',N,'SelEdgeIndex',sH_max_ind,'I',I,'meanY',meanY,'meanG',meanG);
    end
end

% Now apply to entire dataset
switch flags.corrmethod 
    case 1
        IN.TsCovars = G; sY = nk_PartialCorrelationsObj(Y, IN);
    case 2
        sY = nk_PerfNorm(Y, G, meanY, indL);
        
    case 3
        sY = nk_PerfRemMeanDiff(Y, G, G, meanY, meanG);
        
end
if flags.display
    figure; 
    plot(mean(Y,2),'b'); hold on
    plot(mean(sY,2),'r')
end
fprintf('\nDone.\n')