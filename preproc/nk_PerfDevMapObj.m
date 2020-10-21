function [sY, IN, Yrecon] = nk_PerfDevMapObj(Y, IN)

switch IN.DEVMAP.algostr
    case {'pls','spls'}
        devfun = 'PerfPLSObj'; params = IN.DEVMAP;
end

if isfield(params,'mpp') % Apply model to data
    if iscell(Y) 
        sY = cell(1,numel(Y)); Yrecon = sY;
        for i=1:numel(Y), 
            % Define active indices depending on training or testing situation
            if isfield(IN,'TsInd'), iX = params.covmat(IN.TsInd{i},:); else, iX = params.covmat; end
            [sY{i}, ~, Yrecon{i}] = feval ( devfun, iX, Y{i}, params ); 
        end
    else
        if isfield(IN,'TsInd'), X = params.covmat(IN.TsInd,:); else, X = params.covmat; end
        [sY , ~, Yrecon ] = feval( devfun, X, Y, params ); 
    end
else % Train model using data
    if isfield(IN,'TrInd'), 
        X = params.covmat(IN.TrInd,:); 
        if isfield( params,'glabel' ) && ~isempty(params.glabel) 
            X = X(params.glabel(IN.TrInd),:); 
            Y = Y(params.glabel(IN.TrInd),:); 
        end
    else, 
        X = params.covmat; 
        if isfield( params,'glabel' ) && ~isempty(params.glabel)
            X = X(params.glabel,:); 
            Y = Y(params.glabel,:); 
        end
    end
    [ sY, params, Yrecon ] = feval( devfun, X, Y, params); 
end

switch IN.DEVMAP.algostr
    case {'pls','spls'}
        IN.DEVMAP = params;
end

function [sY, IN, Yrecon] = PerfPLSObj(X,Y,IN)

if ~isfield(IN,'mpp') 
    [m,n]=size(Y);
    if isfield(IN,'loo') && IN.loo   
       Yrecon = zeros(m,n); fprintf('\n');
       for i=1:m
           ind = true(m,1); ind(i)=false;
           switch IN.algostr
               case {'pls','spls'}
                   [ ~, ~, ~, INi ] = nk_PLS(Y(ind,:),X(ind,:),[]);
                   [ ~, ~, Yrecon(~ind,:)] = nk_PLS(Y(~ind,:),X(~ind,:),INi);
           end
           fprintf('.')
       end
    end
    switch IN.algostr
        case {'pls','spls'}
            [ ~, ~, ~, IN ] = nk_PLS(Y,X,IN);
            if ~isfield(IN,'loo') || ~IN.loo, [ ~, ~, Yrecon] = nk_PLS(Y,X,IN); end
    end
else
    switch IN.algostr
        case {'pls','spls'}
            [ ~, ~, Yrecon] = nk_PLS(Y,X,IN);      
    end
end

sY = Y - Yrecon;

