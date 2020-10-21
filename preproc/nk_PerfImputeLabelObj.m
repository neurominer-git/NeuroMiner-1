function [sY, IN] = nk_PerfImputeLabelObj(Y, IN)
% =========================================================================
% FORMAT function [sY, IN] = nk_PerfImputeLabelObj(Y, IN)
% =========================================================================
% Performes imputation using either single-subject median replacement, 
% feature-wise mean replacement or multivariate distance-based NN median 
% imputation. If you want to use hamming or euclidean distances you should 
% scale, unit-normalize or standardize the data first, otherwise the
% distance measure will be dominated by high-variance features.
% 
% I/O arguments:
% IN.Y              : feature matrix
% IN.method         : [ singlemean ] single-subject median replacement if
%                                    NaN value in feature block
%                     [ mean ]       NaN value is replaced by mean of value
%                                    within the given feature
% [ dist, dist2, seuclidean, cosine ] NN-based replacement of NaN using
%                     the median of the IN.K most similar instances with
%                     finite values for the given NaN value. Similarity
%                     is determined from IN.X using IN.method. 'seuclidean' 
%                     and 'cosine' require pdist2 which is available through 
%                     the MATLAB statistics toolbox  
% sY                : The imputed label vector
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Nikolaos Koutsouleris, 10/2015
global VERBOSE MODEFL

%% Prepare
if ~isfield(IN,'method') || isempty(IN.method), IN.method = 'dist2'; end
if ~isfield(IN,'k') || isempty(IN.k), IN.k = 7; end
fnan = find(~isfinite(Y));
ffin = isfinite(Y);
S = IN.X(ffin,:);
SL = Y(ffin,:);
sY = Y;
switch IN.method
    case 'seuclidean'
        C = nm_nanstd(S); C(C==0) = min(C(C~=0));
end
switch MODEFL
    case 'classification'
        frepl = 'mode';
    otherwise
        frepl = 'median';
end

if VERBOSE, fprintf('\tImpute missing labels');end
ll=0;            
for i=1:numel(fnan)
   
    Ti = IN.X(fnan(i),:);
        
    % Compute distance metric
    switch IN.method
        case {'manhattan', 'euclidean'}
             D = pdist2(S, Ti, 'euclidean')';
        case 'seuclidean'
            D = pdist2(S, Ti, 'seuclidean', C')';
        case 'mahalanobis'
            %C = nancov(Xj(:,indi_Yi)); C(C==0) = min(C(C~=0));
            D = pdist2(S, Ti, 'mahalanobis')';
        otherwise
            D = pdist2(S, Ti,IN.method)';
    end
    % Sort training cases according to their proximity to the
    % test case whose value will be imputed.
    [Ds, ind] = sort(D,'ascend');
    if size(Ds,1) < IN.k, kx = size(Ds,1); else kx = IN.k; end
    % Here we take the unweighted median of the kx imputation
    % training cases
    mn = feval(frepl, SL(ind(1:kx)));
    sY(fnan(i)) = mn;
    %if VERBOSE,fprintf('+'), end
    ll=ll+1; 
end

if VERBOSE, fprintf('\n\t\t\t%g subject(s) with NaNs labels imputed.',ll); end
