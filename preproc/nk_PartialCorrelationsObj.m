function [sY, IN] = nk_PartialCorrelationsObj( Y, IN )
% =========================================================================
% FORMAT [Y, IN] = nk_PartialCorrelationsObj(Y, IN)
% =========================================================================
% Remove nuisance effects IN.G from Y (opt. using predefined estimators)
%
% I\O Arguments:
% -------------------------------------------------------------------------
% Y                 : M cases x N features data matrix
% IN.G              : The covariate(s) to regress out from Y
% IN.nointercept    : Include an intercept in the model or not
% IN.subgroup       : Index vector of cases in Y to compute the beta(s) from
% IN.beta           : The estimated beta coefficients
% IN.revertflag     : Increase (=1) or attenuate (=0) IN.G effects 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Nikolaos Koutsouleris, 08 / 2020

method = @PartialCorrelationsObj; methodsel = 1;
if exist('IN','var') && ~isempty(IN) && isfield(IN,'METHOD') && IN.METHOD == 2
    method = @CombatObj; methodsel = 2;
end

% =========================== WRAPPER FUNCTION ============================ 
if iscell(Y) && exist('IN','var') && ~isempty(IN)
    sY = cell(1,numel(Y));
    for i=1:numel(Y),
        switch methodsel 
            case 1
                if isfield(IN,'beta') && ~isempty(IN.beta)
                    IN.G = IN.TsCovars{i};
                else
                    IN.G = IN.TrCovars{i};
                end
            case 2
                if isfield(IN,'estimators') && ~isempty(IN.estimators)
                    IN.G = IN.TsCovars{i};
                    IN.M = IN.TsMod{i};
                else
                    IN.G = IN.TrCovars{i};
                    IN.M = IN.TrMod{i};
                end
        end
        [sY{i}, IN] = method (Y{i}, IN); 
    end
else
    switch methodsel 
        case 1
            if isfield(IN,'beta') && ~isempty(IN.beta)
                if iscell(IN.TsCovars)
                    IN.G = IN.TrCovars;
                else
                    IN.G = IN.TsCovars;
                end
            else
                IN.G = IN.TrCovars;
            end

        case 2
            if isfield(IN,'estimators') && ~isempty(IN.estimators)
                if iscell(IN.TsMod)
                    IN.G = IN.TrCovars;
                    IN.M = IN.TrMod;
                else
                    IN.G = IN.TsCovars;
                    IN.M = IN.TsMod;
                end
            else
                IN.G = IN.TrCovars;
                IN.M = IN.TrMod;
            end
    end

    [ sY, IN ] = method( Y, IN );
end

% =========================================================================
function [Y, IN] = PartialCorrelationsObj( Y, IN )

if isempty(IN),eIN=true; else, eIN=false; end

if eIN|| ~isfield(IN,'G') || isempty(IN.G), error('No covariates defined in parameter structure'), end

if eIN || (~isfield(IN,'nointercept') || isempty(IN.nointercept) || ~IN.nointercept ) 
     %Create intercept vecotr
    interceptflag = true;
    intercept = ones(size(IN.G,1),1);
    %Check if intercept is already included in matrix to avoid double
    %intercept removal
    if isfield(IN,'beta') && ~isempty(IN.beta)
        if size(IN.beta,1) == size(IN.G,2), interceptflag = false; end
    end
else
    interceptflag = false;
end

if interceptflag
    %fprintf(' ... adding intercept to covariate matrix')
    IN.G = [intercept IN.G];
end

if eIN || ~isfield(IN,'beta') || isempty(IN.beta), 
    if ~isfield(IN,'subgroup') || isempty(IN.subgroup)
        % Compute IN.beta from entire population
        IN.beta = pinv(IN.G) * Y; 
    else
        % Compute IN.beta from a subgroup of observations
        IN.beta = pinv(IN.G(IN.subgroup,:)) * Y(IN.subgroup,:);
    end
end
if eIN || ~isfield(IN,'revertflag') || isempty(IN.revertflag) || ~IN.revertflag
    Y = Y - IN.G * IN.beta;
else
    Y = Y + IN.G * IN.beta;
end

% =========================================================================
function [ Y, IN ] = CombatObj( Y, IN )

if isempty(IN),eIN=true; else, eIN=false; end

if eIN|| ~isfield(IN,'G') || isempty(IN.G), error('No covariates defined in parameter structure'), end

% Combat needs data, batch variables and covariates transposed 
Y=Y'; G = IN.G'; if islogical(G), [~,G] = max(G); end
M = []; if isfield(IN,'M'), M = IN.M; end

if eIN || (~isfield(IN,'estimators') || isempty(IN.estimators) ) 
    
    [~,IN.estimators] = combat( Y, G, M );
end

Y = combat( Y, G, M, IN.estimators);

Y=Y';