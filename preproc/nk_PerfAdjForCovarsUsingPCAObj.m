function [sY, IN, dT] = nk_PerfAdjForCovarsUsingPCAObj(Y, IN, S)
% =========================================================================
% function [adjT, IN] = nk_AdjForCovarsUsingPCAObj(T, IN)
% =========================================================================
% ----- INPUTS -----
% IN.S :            Source Matrix (row = patterns, cols = features)
%                   Used to (1) determine which eigenvectors are correlated
%                   with CS, and (2) to reduce and reconstruct T
% IN.G :            Covariate matrix 
%
% Optional fields in IN include:
% IN.recon :        back-project adjusted target matrix to input space 
%                   (def: true)
% IN.varop :        operator for identification of correlated factors 
%                   (def: gt)
% IN.corrmeth :     identification method (1 = pearson (def), 
%                   2 = spearman, 3 = anova)
% IN.corrthresh :   cutoff for identification (def: 0.3)
% IN.DR :           Dimensionality reduction parameter structure (see
%                   nk_DimRed_config)
% IN.indX :         logical index vector to subcohort of IN.S
%
% Y :               Target matrix, from which the covariate effects will be 
%                   removed
%
% ----- OUTPUT -----
% adjT :            Adjusted target matrix 
% IN :              Parameter structure containing orig. and comp. params
% =========================================================================
% (c) Nikolaos Koutsouleris, 1/2016

% =========================== WRAPPER FUNCTION ============================ 
if ~exist('S','var'), S=[]; end
if iscell(Y) && exist('IN','var') && ~isempty(IN)
    sY = cell(1,numel(Y)); dT = [];
    for i=1:numel(Y), [sY{i}, IN] = PerfAdjForCovarsUsingPCAObj(Y{i}, IN, S); end
else
    [ sY, IN, dT ] = PerfAdjForCovarsUsingPCAObj( Y, IN, S );
end

end

% =========================================================================
function [adjT, IN, dT] = PerfAdjForCovarsUsingPCAObj(T, IN, S)
global VERBOSE

% Check existence of paramater structure
if ~exist('IN','var') || isempty(IN),             
    error('No parameter structure specified! Abort!');  
end
if ~isfield(IN,'recon') || isempty(IN.recon),
    IN.recon = true;
end
if ~isfield(IN,'varop') || isempty(IN.varop)
    IN.varop = 'gt';
end
    
compfl = false;
% Check whether variance removal parameters have been already computed
if ~isfield(IN,'ind0') || isempty(IN.ind0)
    
    compfl = true;
    % The source data matrix to compute the factorized matrix from
    if ~exist('S','var') || isempty(S), 
        error('No training matrix specified in parameter structure'); 
    end
    % The source covariate matrix to compute correlation coefficients with
    % factorized data matrix
    if (~isfield(IN,'G') || isempty(IN.G)),
        error('No target vector / matrix specified in parameter structure'); 
    end
    % The defaults correlation method
    if ~isfield(IN,'corrmeth') || isempty(IN.corrmeth),
        IN.corrmeth = 1;                                    
    end
    switch IN.corrmeth
        case 1
            corrmeth = 'pearson';
        case 2
            corrmeth = 'spearman';
        case 3
            corrmeth = 'anova';
    end

    % The default correlation strength 
    if ~isfield(IN,'corrthresh') || isempty(IN.corrthresh),
        IN.corrthresh = 0.3;                                
    end
    
    % if not otherwise specified use the entire source data matrix
    if ~isfield(IN,'indX') || isempty(IN.indX)
        IN.indX = true(size(S,1),1);
    end
    
    % if not otherwise specified use PCA in the 'percentage of
    % sum(eigenvalues) mode'
    if ~isfield(IN,'DR') || isempty(IN.DR) || ~isfield(IN.DR,'RedMode')
        IN.DR.DRsoft = 1; 
        IN.DR.RedMode = 'PCA';
        IN.DR.PercMode = IN.dimmode;  
    end
    
    % Run dimensionality reduction on the source matrix
    [dS, IN] = nk_PerfRedObj(S(IN.indX,:),IN);
    IN.C = zeros(size(dS,2),size(IN.G,2));
    
    % Identify correlated factors in the factorized source matrix
    switch corrmeth
        case {'pearson', 'spearman'}
            for i = 1:size(IN.G,2)
                % Determine correlations with covars in the source matrix
                if VERBOSE, fprintf('\nWorking on covariate #%g', i); end
                IN.C(:,i) = abs(nk_CorrMat(dS,IN.G(IN.indX,i),corrmeth)');
            end
        case 'anova'
            % This is the more powerful option if we have multiple covars
            % for which we would like to identify respective variance
            % components
            RES.X = [ones(size(IN.G(IN.indX,:),1),1) IN.G(IN.indX,:)];
            RES = nk_PerfANOVAObj(dS, RES);
            IN.C = sqrt(RES.R2);
    end

    % Threshold eigenvariate correlations with covars
    IN.subthresh = single(feval(IN.varop, IN.C, IN.corrthresh));
    
    % Check whether correlated factors exist or not!
    if ~sum(IN.subthresh),
        IN.ind0 = true(1,size(dS,2));
        warning('\nNo variance components identified at %g threshold. Returning unchanged matrix!', IN.corrthresh);
        switch IN.recon
            case 1
              adjT = T;
            case 2
              if isempty(T) ||(isequal(T, S) && sum(IN.indX) == size(S,1))
                 adjT = dS; 
              else
                 adjT = nk_PerfRedObj(T,IN);
              end
        end
        dT = adjT;
        return
    end

    % Determine which eigenvariates have to be extracted
    IN.ind0 = ~sum(IN.subthresh,2);
end
if ~any(IN.ind0) 
    warning('All eigenvariates meet threshold criterion [%s %g]!\nCheck your data and your settings.',IN.varop,IN.corrthresh);
end
    
% Now project target matrix to source matrix PCA space
if VERBOSE, fprintf('\nProjecting target matrix to source PCA space'); end
if isempty(T) ||(exist('S','var') && isequal(T, S) && sum(IN.indX) == size(S,1))
    if exist('dS','var')
        dT = dS; 
    else
        % Run dimensionality reduction on the source matrix
        dT = nk_PerfRedObj(S(IN.indX,:),IN);
    end
else
    dT = nk_PerfRedObj(T,IN);
end

switch IN.recon
    case 1
        % Reconstruct target matrix only with / without extracted variance
        % components
        if VERBOSE, fprintf('\nReconstructing target matrix based on selected components.'); end
        tT = bsxfun(@plus, IN.mpp.vec(:,IN.ind0)* dT(:,IN.ind0)' , IN.mpp.sampleMean')';
        adjT = zeros(size(dT,1),numel(IN.indNonRem)); adjT(:,IN.indNonRem) = tT;
    case 2
        if VERBOSE, fprintf('\nLimiting target matrix to selected components.'); end
        % Limited projected target matrix to selected variance components
        adjT = dT(:,IN.ind0);
end

if VERBOSE, 
    if compfl
        fprintf(['\nProcessing finished. %g/%g components with r %s %g (max: %g) were identified.', ...
            '\nRespective variance components was extracted from the input matrix.\n'],sum(~IN.ind0),size(IN.mpp.vec,2), IN.varop, IN.corrthresh, max(IN.C(~IN.ind0))); 
    else
        fprintf('\nProcessing finished. %g variance components were removed from input matrix', sum(~IN.ind0))
    end
end

end