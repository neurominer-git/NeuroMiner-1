function [ts, rs, ds, Model, map] = nksp_GetTestPerf(Xtest, Ytest, Features, Model, X, nonevalflag)
% =====================================================================================
% function [ts, rs, ds] = nk_GetTestPerf(Xtest, Ytest, Features, Model, X, nonevalflag)
% =====================================================================================
% Applies Model to Xtest and optionally (nonevalflag) evaluates prediction
% performance using labels of the test data. Training data have to be
% provided for algorithms that need them for evaluating the test data
% Inputs:
% Xtest :       The test data (row=obs, cols=feats)
% Ytest :       The test labels
% Features :    A feature mask contained which can be either (a) a numeric 
%               array with 0 and non-0 entries indicating features not to be 
%               used for prediction, or (b) a cell array with 
%               numel(F) = numel(Model). Each cell is is feature mask array 
%               as defined in (a).
% Model :       A model container with either (a) a single model structure, 
%               or (b) a cell array of models structures.
%
% Output: 
% ts :          Test performance
% rs :          Predicted labels: 
%                   in classification:  rs = sign(ds),
%                   in regression:      rs = ds;
% ds :          Prediction scores
% =====================================================================================
% (c) Nikolaos Koutsouleris, 12/2019

global PREDICTFUNC EVALFUNC SVM MODEFL TRANSPARENT

% ******************************* Prepare *********************************
s=1;
if ~isempty(Features)
    if iscell(Features), 
        s=size(Features{1},2); 
    else
        s=size(Features,2); 
    end
end
ts = zeros(s,1);
if iscell(Xtest),
    % Check and remove Nan cases
    [Xtest{1}, Ytest, I] = nk_ManageNanCases(Xtest{1}, Ytest);
    nsubj = size(Xtest{1},1); 
else, 
    [Xtest, Ytest, I] = nk_ManageNanCases(Xtest, Ytest);
    nsubj = size(Xtest,1); 
end
rs = zeros(nsubj,s); ds = rs;
if ~exist('nonevalflag','var'), nonevalflag = false; end

for k=1:s % Loop through all feature subspaces
        
    % ************** Get test data using current feature mask *************
    tXtest  = nk_ExtractFeatures(Xtest, Features, [], k);
    
    % ******************** Get Model of current subpace *******************
    if iscell(Model), md = Model{k}; else md = Model; end
    
    % Check test matrix for NaN observations and remove them
    switch SVM.prog
        case 'SEQOPT'
            [ rs(:,k), ds(:,k), md ] = feval(PREDICTFUNC, tXtest, Ytest, md);
            if iscell(Model), Model{k} = md; else, Model = md; end
            if ~nonevalflag, ts(k) = feval(EVALFUNC, Ytest, ds(:,k)); end
        case 'WBLCOX'
            % we will treat the Weibull-Cox regression model differently
            [ rs(:,k), ds(:,k), md ] = nk_GetTestPerf_WBLCOX(tXtest, md);
            if iscell(Model), Model{k} = md; else, Model = md; end
            if ~nonevalflag, ts(k) = feval(EVALFUNC, Ytest, rs(:,k)); end
        otherwise
            [ rs(:,k), ds(:,k) ] = feval(PREDICTFUNC, X, tXtest, Ytest, md, Features, k);
            % Adjust probabilities if probabilistic output has been
            % geenerated by the probabilistic algorithm
            if SVM.RVMflag && ~strcmp(MODEFL,'regression') , ds(:,k) = nk_CalibrateProbabilities(ds(:,k)); end
            % Return performance measure as defined by EVALFUNC
            if ~nonevalflag, ts(k) = feval(EVALFUNC, Ytest, ds(:,k)); end
    end
    
    %% LIME EVALUATION
    if TRANSPARENT.flag
        switch TRANSPARENT.algo
            case 'lime'
                map = sp_Lime()
            case 'sharp'
        end
    end
    

end

% Check and add-back Nan cases
[rs, ds] = nk_ManageNanCases(rs, ds, I);