function [sY, IN] = nk_PerfPLSObj(X, Y, IN)

if iscell(Y) && exist('IN','var') && ~isempty(IN)
    sY = cell(1,numel(Y)); 
    for i=1:numel(Y), 
        % Define active indices depending on training or testing situation
        if isfield(IN,'TsInd'), IN.indY = IN.TsInd{i}; else IN.indY = []; end
        sY{i} =  PerfNormObj(Y{i}, IN ); 
    end
else
    % Define active indices depending on training or testing situation
    if isfield(IN,'meanY') 
       if isfield(IN,'TsInd'), IN.indY = IN.TsInd; else IN.indY = []; end
    else
       if isfield(IN,'TrInd'), IN.indY = IN.TrInd; else IN.indY = []; end
    end
    [ sY, IN ] = PerfPLSObj( Y, IN );
end

function [Yrecon, IN] =  PerfPLSObj(X,Y,IN)

if ~isfield(IN,'BETA')
    [m,n]=size(Y);
    if isfield(IN,'loo') && IN.loo   
       Yrecon = zeros(m,n); fprintf('\n');
       for i=1:m
           ind = true(m,1); ind(i)=false;
           [~,~,~,~,IN.BETA] = plsregress(X(ind,:),Y(ind,:),IN.ncomp);
           Yrecon(~ind,:) = [1 X(~ind,:)]*IN.BETA;
           fprintf('.')
       end
    end
    [IN.XL,IN.YL,IN.XS,IN.YS,IN.BETA,IN.PCTVAR] = plsregress(X,Y,IN.ncomp);
else
    m=size(X,1);
    Yrecon = [ones(m,1) X]*IN.BETA;
end

