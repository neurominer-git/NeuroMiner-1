% =========================================================================
% FORMAT function [param, model] = nk_GetParam_LIBSVM(Y, label, ModelOnly, 
%                                                                 ...cmdstr)
% =========================================================================
% Train LIBSVM models and evaluate their performance using Y & label, 
% SlackParam, KernParam
% if ModelOnly = 1, return only model
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Nikolaos Koutsouleris, 07/2015

function [param, model] = nk_GetParam2_CCSSVM(Y, label, ModelOnly, cmd)
                                            
global SVM EVALFUNC LIBSVMTRAIN MODEFL CMDSTR                         

param = [];flw = 0;
if strcmp(LIBSVMTRAIN,'svmtrain312'), flw = 1; end

% Check if weighting is necessary
cmd = nk_SetWeightStr(SVM, MODEFL, CMDSTR, label, cmd); cmd = [cmd CMDSTR.quiet];

% Check if sample weighting is necessary (currently regression only)
if flw
    W = ones(numel(label),1);
    if strcmp(MODEFL,'regression') && SVM.LIBSVM.Weighting
        W = nk_WeigthDataInstanceHisto(label);
    end
end

if iscell(Y) 
   
    % MKL-based learning not implemented yet
   
else % Univariate case
    
    if size(label,1) ~= size(Y,1), label = label'; end
   
    if isfield(SVM,'AdaBoost') && SVM.AdaBoost.flag && flw
        
        N = length(label); % X training labels
        W = 1/N * ones(N,1); %Weights initialization
        for m=1:SVM.AdaBoost.BoostIter
            %Calculate the error and alpha in adaBoost with cross validation
            model               = svmtrain312( W, label, Y, cmd ); 
            Xout                = svmpredict312( label, Y, model );
            s1                  = sum( (label ==  1) .* (Xout) .* W);
            s2                  = sum( (label == -1)  .* (Xout) .* W);
            if s1 == 0 && s2 == 0, break; end
            % Compute Alpha
            alpha              = 0.5*log((s1 + eps) / (s2+eps));  
            % update the weight
            W                   = exp( -1 * (label .* Xout .* alpha));
            W                   = W/norm(W);
        end  
    else
        if flw
            model = feval( LIBSVMTRAIN, W, label, Y, cmd);
        else
            model = feval( LIBSVMTRAIN, label, Y, cmd );
        end
    end
    if ~ModelOnly
        [param.target, param.dec_values] = nk_GetTestPerf_LIBSVM([], Y, label, model) ;
        if SVM.RVMflag, param.dec_values = nk_CalibrateProbabilities(param.dec_values); end
        param.val = feval(EVALFUNC, label, param.target);
    end

end