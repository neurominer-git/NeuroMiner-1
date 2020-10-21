function [pY, mpp] = nk_PerfICA(Y, mpp, k)

% Default values for 'pcamat' parameters
Dim               = size(Y,1);
firstEig          = 1;
lastEig           = Dim;
interactivePCA    = 'off';
verbose           = 'on';

if ~exist('mpp','var') || isempty(mpp)

    %  Calculate PCA
    options.ReducedDim = 0;
    [mpp.eigvector, mpp.eigvalue, mpp.sampleMean] = PCA(Y, options);
    Y = bsxfun(@minus, Y, mpp.sampleMean) * mpp.eigvector;
    
    % Whiten data
    [mpp.whitesig, mpp.whiteningMatrix, mpp.dewhiteningMatrix] = whitenv ...
						     (Y, mpp.E, mpp.D, verbose);
                         
    % Apply ICA
    [mpp.IC, mpp.A, mpp.W] = fastica(Y, 'approach', 'symm', 'NumOfIC', k, ...
                                        'pcaE', mpp.E, ...
                                        'pcaD', mpp.D, ...
                                        'whiteSig', mpp.whitesig, ...
                                        'whiteMat', mpp.whiteningMatrix, ...
                                        'dewhiteMat', mpp.dewhiteningMatrix);
else
    % Remove mean in out-of-sample data
    Y = bsxfun(@minus, Y, mpp.SampleMean);
    % Whiten data
    [mpp.whitesig, mpp.whiteningMatrix, mpp.dewhiteningMatrix] = whitenv ...
						     (Y, mpp.E, mpp.D, verbose);
    % Apply ICA
    [~, mpp.W] = fastica(Y, 'approach', 'symm', 'NumOfIC', k, ...
                                        'pcaE', mpp.E, ...
                                        'pcaD', mpp.D, ...
                                        'whiteSig', mpp.whitesig, ...
                                        'whiteMat', mpp.whiteningMatrix, ...
                                        'dewhiteMat', mpp.dewhiteningMatrix);
end

pY = mpp.W';