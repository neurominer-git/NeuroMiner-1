function [adjT, IN] = nk_AdjForCovarsUsingPCA(T, IN, S, CS, corrthresh)
% =========================================================================
% function [adjT, Param] = nk_AdjForCovarsUsingPCA(T, IN, S, CS, corrthresh)
% -------------------------------------------------------------------------
% Inputs:
% S :               Source Matrix (row = patterns, cols = features)
%                   Used to (1) determine which eigenvectors are correlated
%                   with CS, and (2) to reduce and reconstruct T
% CS :              Covariate matrix 
% T :               Target matrix, from which the covariate effects will be 
%                   removed
% Output:
% adjT :
% IN :
% =========================================================================
% (c) Nikolaos Koutsouleris, 12/2015
global VERBOSE

adjT = zeros(size(T));

if ~exist('IN','var') || isempty(IN)  
    
    if ~exist('corrthresh','var') || isempty(corrthresh)
        IN.corrthresh = 0.3;
    else
        IN.corrthresh = corrthresh;
    end
    IN.DR.DRsoft= 1;
    IN.DR.RedMode= 'PCA';
    IN.DR.PercMode = 3;
    IN.DR.dims= 0.9;
    % Run PCA on the source matrix
    [dS, IN] = nk_PerfRedObj(S,IN);
    IN.C = zeros(IN.DR.dims,size(CS,2));

    for i = 1:size(CS,2)
        % Determine correlations with covars in the source matrix
        if VERBOSE, fprintf('\nWorking on covariate #%g', i); end
        IN.C(:,i) = abs(nk_CorrMat(dS,CS(:,i),'spearman'));
    end

    % Threshold eigenvariate correlations with covars
    IN.supthresh = IN.C < IN.corrthresh;

    % Determine eigenvariate for removal
    IN.ind0 = logical(prod(IN.supthresh,2));
end

% Now project target matrix to source matrix PCA space
if VERBOSE, fprintf('\nProjecting target matrix to source PCA space'); end
dT = nk_PerfRedObj(T,IN);

% Reconstruct target matrix without eigenvariates
if VERBOSE, fprintf('\nReconstructing target matrix without identified variance components.'); end
tT = bsxfun(@plus, IN.mpp.vec(:,IN.ind0)* dT(:,IN.ind0)' , IN.mpp.sampleMean')';

adjT(:,IN.indNonRem) = tT;

if VERBOSE, fprintf(['\nProcessing finished. %g/%g eigenvariates with correlations >%g were identified.', ...
    '\nRespective variance was removed from the data.\n'],sum(~IN.ind0),IN.DR.dims, IN.corrthresh); end

end