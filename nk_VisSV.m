function [Vmean, V] = nk_VisSV(model, X, Y)
% 
% function [Vmean, V] = nk_VisSV(model, X, Y)
% 
% Computes the minimum difference vector matrix V for a given SVM model 
% (currently only LIBSVM) by finding the minimum difference vector for each SV 
% of class +1 to the SVs of class -1.
%
% INPUTS
% * model:        LIBSVM model (currently)
% * X:            Preprocessed data (MxN matrix, m = obs., n = features
% * Y:            Labels or Target values
%
% OUTPUTS
% * V:            Matrix with min. diff. vector for each SV of class +1
%                 to the SVs of class -1
% * Vmean:        Mean of V across minimum SV difference vectors of class +1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Nikolaos Koutsouleris, 12/2010


% Check whether we have a regression or classification model

if length(model.nSV) > 1
    
    %%%%%%%%%%%%%%%%%%%%%%% CLASSIFICATION MODEL %%%%%%%%%%%%%%%%%%%%%%%%%%
	% nSV1,2 = number of SV of class 1,2
    SVind     = ismember(X, model.SVs, 'rows'); 
    indSV1  = ( Y == 1    &   SVind == 1 ); % index to SVs of class +1
    indSV2  = ( Y == -1   &   SVind == 1 ); % index to SVs of class -1
    nSV1    = sum(indSV1);	nSV2 = sum(indSV2);
    oSV = nSV1; iSV = nSV2; oInd = indSV1; iInd = indSV2; %flg=1;
    
    % Initialize
    V = zeros(oSV,size(X,2));
    
	% 1) Calculate SV difference vectors for each SV1_j and SV2_k=1...n
    SV1 = X(oInd,:); % Support vectors of outer loop (class +1) 
    SV2 = X(iInd,:); % Support vectors of inner loop (class -1)
    
    for j=1:oSV % Loop through the SVs of class +1
        minSVmeandiff = Inf;
        for ki=1:iSV % Loop through the SVs of class -1
            % Search for the min. diff. vector
			% SVdiff                  = (SV1(j,:) - SV2(k,:))*flg;
            % SVmeandiff              = abs(mean(SVdiff));
            SVdiff     = (SV1(j,:) - SV2(ki,:)); SVmeandiff = sum(abs(SVdiff));
            if SVmeandiff < minSVmeandiff, V(j,:) = SVdiff; minSVmeandiff = SVmeandiff; end
        end
    end
    
else
    %%%%%%%%%%%%%%%%%%%% REGRESSION / ONE-CLASS MODEL %%%%%%%%%%%%%%%%%%%%%
    % 1) Simply compute the average SV vector across all SVs in the model
    indSV       = ismember(X, model.SVs, 'rows');
    V           = mean(X(indSV,:));
   
end

% 2) Compute the mean difference vector across D (dimensions)
if size(V,1)>1,Vmean = mean(V); else Vmean = V; end

return



