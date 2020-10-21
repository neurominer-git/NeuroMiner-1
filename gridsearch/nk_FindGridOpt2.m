function [ipos, ind0] = nk_FindGridOpt2(tr, ts, c, action, perc)
%
% function [xpos, ypos] = nk_FindGridOpt(tr, ts, c, action, alttest)
% 
% This function determines the optimum C and kernel parameters after the
% grid computations have been finished. It can be used for binary /
% multi-class or regression model optimization
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Nikolaos Koutsouleris 11/2015

global GRD SVM VERBOSE

otr = tr; ots = ts;

switch action
    case 'min'
        fmult = -1; evalop = 'le';
    case 'max'
        fmult= 1; evalop = 'ge';
end
% Apply regularization parameters to grid in order to improve
% generalization performance. We look for optima with low model complexity
% (or alternatively high ensemble diversity).
% This is achived by adding penality to optima with high complexity (or low
% diversity) using the functions: 
%
% Model complexity:     opt = opt - lambda*(n_sv/n_obs).^Gamma;
%   , where n_sv = Number of support (relevance) vectors and n_obs = Number
%     of samples in the training population
%
% Ensemble diversity:   opt = opt + lambda*

if isfield(GRD,'OptRegul') && GRD.OptRegul.flag
    
    fac = fmult*GRD.OptRegul.lambda*(c.^GRD.OptRegul.big_gamma);
    
    switch SVM.GridParam
        case {5, 6, 7, 13, 15}
            fac = fac./100;
    end
    tr = tr - fac; 
end

ind0 = true(size(tr));
imx = feval(action,tr);
[mx,im] = feval(action,tr(ind0));

if isfield(GRD,'NodeSelect')

    switch GRD.NodeSelect.mode

        case 1 % Single-Node selection
            ind0 = ind0 & tr == mx;
            if sum(ind0)>1 
                if (exist('c','var') && ~isempty(c)) && numel(unique(c)) > 1
                    mxc = min(c(ind0));
                    indc = c == mxc;
                    ind0 = ind0 & indc;
                    if sum(ind0) > 1, 
                        elimind = find(ind0); % take only first
                        ind0(elimind(2:end)) = false;
                    end
                else
                    elimind = find(tr==mx); % take only first
                    ind0(elimind(2:end)) = false;
                end
            end

        case {2, 4} % Multi-Node selection (predetermined threshold)

            if GRD.NodeSelect.mode == 2
                percX = GRD.NodeSelect.perc;
            else
                percX = perc;
            end
            thresh = percentile(tr(ind0), percX);
            ind0 = ind0 & feval(evalop, tr,thresh);
            fprintf('\nSelected %g parameter nodes above the %g%% percentile.', sum(ind0(:)),percX)
    end
end

ipos = find(ind0 == true);

if VERBOSE;
    if isfield(GRD,'OptRegul') && GRD.OptRegul.flag   
        if isfield(GRD.OptRegul,'type')
            switch GRD.OptRegul.type
                case 1
                    strtype = ['Sparseness [' GRD.OptRegul.RegulTypeComplexity ']'];
                case 2 
                    strtype = ['Diversity [' GRD.OptRegul.RegulTypeDiversity ']'];
                case 3
                    strtype = ['Average of Sparseness [' GRD.OptRegul.RegulTypeComplexity '] & Diversity [' GRD.OptRegul.RegulTypeDiversity ']' ];
            end
        else
            strtype = 'Sparseness';
        end
        fprintf('\nRegularization (%s) of model selection enabled: Lambda = %g, Big-Gamma = %g', ...
            strtype, GRD.OptRegul.lambda, GRD.OptRegul.big_gamma);
        fprintf('\nSearching for CV1-%s that meets criteria.', action);
        if ~isempty(mx)
            fprintf('\nInitial CV1-%s score = %1.2f',action,otr(im))
            fprintf(' (regularization score = %1.2f)',mx)
        else
            fprintf('\nInitial CV1-%s score = %1.2f',action,imx)
        end
        fprintf('\nOptimum at CV1 = %1.2f, CV2 = %1.2f',otr(ipos(1)),ots(ipos(1)))
        if  (exist('c','var') && ~isempty(c)),fprintf('\nComplexity at optimum: %1.2f%%',c(ipos(1))); end
    else
        fprintf('\nRegularization of model selection disabled.'); 
    end
else
    if isfield(GRD,'OptRegul') && GRD.OptRegul.flag, fprintf('\nRegularization of model selection enabled.'); end 
end
