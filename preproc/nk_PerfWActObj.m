function [ sY, IN ] = nk_PerfWActObj( Y, IN )
% =========================================================================
% FORMAT function [sY, IN] = nk_PerfWActObj(Y, IN)
% =========================================================================
% This function applies precomputed weight vector(s) over the features to
% the predictor data either by multiplying features with their weights (soft
% feature selection with a given exponent) or by selecting features 
% above given percentile threshold (hard feature selection). 
% In the latter case, the suprathreshold data can also be clustered if it 
% is in 3D format (e.g. neuroimaging data).
% 
% Inputs/Outputs: 
% -------------------------------------------------------------------------
% Y                   : M cases x N features data matrix
% IN                  : Input parameter structure 
%   W_ACT.clustflag   : tells the script to cluster the data   
%   W_ACT.opt         : contains the current hyperparameter for the script
%                       which can be exponents for soft feature selection 
%                       or percentile thresholds for hard feature selection
%   W_ACT.threshvec   : contains all percentiles defined by the user
%   Thresh            : Hard feature selection threshold computed using the
%                       based on W_ACT.opt
%   Mask              : A pre-defined mask for feature selection where each
%                       value in the mask will be used to summarize the data 
%                       in the respective label (see nk_ExtractClusterData.m). 
%                       Instead of using a pre-defined mask the 
%                       script can call the clustering routine of SPM to 
%                       compute a mask based on the thresholded data.
%   W                 : The n*p weight vector/matrix defining feature 
%                       relevance, where n is the feature dimensionality.
%                       If p>1 then the feature selection operation will be
%                       performed n times with each weight vector in the
%                       matrix and the results will be concatenated into sY
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Nikolaos Koutsouleris, 08/2020

% =========================== WRAPPER FUNCTION ============================ 
if iscell(Y) && exist('IN','var') && ~isempty(IN)
    sY = cell(1,numel(Y)); 
    for i=1:numel(Y), [sY{i}, IN] = PerfWActObj(Y{i}, IN); end
else
    [ sY, IN ] = PerfWActObj(Y, IN );
end

% =========================================================================

function [ nY, IN ] = PerfWActObj( Y, IN )
global VERBOSE

% Defaults
if isempty(IN),          error('Suitable input structure is missing. See the functions'' help for more information.'); end
if ~isfield(IN,'W_ACT'), error('Weighting parameter structure is missing! Add a W_ACT substructure to your input structure.'); end
if ~isfield(IN,'W'),     error('Weighting vector is missing! Add a weighting vector W to your input structure.'); end
if ~isfield(IN,'Mask'), IN.Mask = []; end

nW = size(IN.W,2);
nY = [];

%Loop through multiple weight vectors and expand feature space if needed
for i=1:nW 
    
    if isfield(IN.W_ACT,'opt')
        Params_desc = IN.W_ACT.Params_desc;
        opt = IN.W_ACT.opt;
    else
        error('\nNo percentile threshold has been provided to extract futures.')
    end

    switch IN.W_ACT.softflag 
        
        % Soft feature selection by multiplying features with respective
        % weights
        case 1 
            % Extract parameters for exponential multiplier of W
            if ~isfield(IN,'Weights') || isempty(IN.Weights)
                if i==1
                    IN.Weights = zeros(size(IN.W));
                    IN.ind = true(size(IN.Weights));
                end
                ExpMult = nk_ReturnParam('ExpMult',Params_desc, opt); 
                IN.Weights(:,i) = IN.W(:,i).^ExpMult; 
            end
            % Soft feature selection
            Y = Y .* IN.Weights(:,i)';
            if VERBOSE, fprintf('\tWeighting F'); end
        
        % Hard feature selection by selecting feature above the given
        % weight threshold
        case 2
            % extract thresholds
            if ~isfield(IN,'Thresh') || isempty(IN.Thresh)
                if i==1, IN.ind = false(size(IN.W)); end
                % Check whether IN.Thresh exists.
                % If not extract percentile(s) from optimization structure exists
                t = nk_ReturnParam('Thresholds',Params_desc, opt); 
                IN.Thresh(i) = percentile(IN.W(:,i), t);
            end

            if IN.W_ACT.clustflag == 1
                % Here we clusterize suprathreshold features if the input
                % data is in 3D format
                if isempty(IN.Mask)
                    Wthresh = zeros(size(IN.W(:,i)));
                    IN.ind(:,i) = IN.W(:,i) > IN.Thresh(i);
                    Wthresh(IN.ind(:,i)) = IN.W(IN.ind(:,i));
                    IN.WMask{i} = nk_Cluster(Wthresh, IN.W_ACT);
                else
                    IN.WMask{i} = IN.Mask;
                end
                Y = nk_ExtractClusterData(Y, IN.WMask{i});
                if VERBOSE, fprintf('\tClusterizing F into %g mean cluster values', nk_Range(IN.WMask{i})); end
            else
                % Here we simply extract features above the given threshold
                IN.ind(:,i) = IN.W(:,i) >= IN.Thresh(i);
                Y = Y(:, IN.ind(:,i)); 
                if VERBOSE, fprintf('\tSelecting %g / %g feats at %g', size(Y,2), numel(IN.W(:,i)), IN.Thresh(i)); end
            end

    end
    nY = [nY Y];
end
