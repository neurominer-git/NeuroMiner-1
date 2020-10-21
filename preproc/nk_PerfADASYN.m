function [ YY, LL, CC, fl ] = nk_PerfADASYN(Y, L, IN, C, noconcatfl)

adasyn_beta                     = [];   %let ADASYN choose default
adasyn_kDensity                 = [];   %let ADASYN choose default
adasyn_kSMOTE                   = [];   %let ADASYN choose default
adasyn_normalized               = false;%false lets ADASYN handle normalization
cfl                             = false;
fl                              = true;
L(L==-1)=0;

if exist('IN','var') && ~isempty('IN')
    if isfield(IN,'beta') && ~isempty(IN.beta), adasyn_beta = IN.beta; end
    if isfield(IN,'kDensity') && ~isempty(IN.kDensity), adasyn_kDensity = IN.kDensity; end
    if isfield(IN,'kSMOTE') && ~isempty(IN.kSMOTE), adasyn_kSMOTE = IN.kSMOTE; end
    if isfield(IN,'normalized') && ~isempty(IN.normalized), adasyn_normalized = IN.normalized; end
end

if exist('C','var') && ~isempty(C)
    cfl = true; nC = size(C,2); 
else
    CC = []; 
end

if ~exist('noconcatfl','var') || isempty(noconcatfl)
    noconcatfl = false;
end

L(L==2)=0;
    
if iscell(Y)
    
    N = numel(Y);
    LL = cell(N,1);
    YY = cell(N,1);
    if cfl, CC = cell(N,1); end
    
    % Loop through training data shelves
    for i=1:N
        
        % If covariates are available add them to training matrix
        if cfl,
            nY = size(Y{i},2); tY = [Y{i} C];
        else
            tY = Y{i};
        end
        
        % Perform ADASYN
        [ tYsyn, Lsyn ] = ADASYN(tY, L, adasyn_beta, adasyn_kDensity, adasyn_kSMOTE, adasyn_normalized); 
        
        % Check whether covariates were integrated in training data and
        % split synthetic data into synthetic training and covariate matrices
        if cfl 
            Ysyn = tYsyn(:,1:nY); Csyn = tYsyn(:,nY+1:nY+nC);
            if ~noconcatfl
                CC{i} = [C; Csyn];
            else
                CC{i} = Csyn;
            end
        else
            Ysyn = tYsyn;
        end
        
        % Add synthetic training data to original training data
        if ~noconcatfl,
            YY{i} = [Y{i}; Ysyn]; LL{i}=[L; Lsyn]; 
        else
            YY{i} = Ysyn; LL{i}= double(Lsyn); 
        end
        LL{i}(~LL{i})=-1;

    end
else
    
    % If covariates are available add them to training matrix
    if cfl,
        nY = size(Y,2); tY = [Y C];
    else
        tY = Y;
    end
    
    % Perform ADASYN
    [tYsyn, Lsyn] = ADASYN(tY, L, adasyn_beta, adasyn_kDensity, adasyn_kSMOTE, adasyn_normalized);
    
    % Check whether covariates were integrated in training data and
    % split synthetic data into synthetic training and covariate matrices
    if ~isempty(tYsyn)
        if cfl 
            Ysyn = tYsyn(:,1:nY); Csyn = tYsyn(:,nY+1:nY+nC);
            if ~noconcatfl,
                CC  = [C; Csyn];
            else
                CC = Csyn;
            end
        else
            Ysyn = tYsyn;
        end
        if ~noconcatfl
            YY = [Y; Ysyn]; LL = [L; Lsyn]; 
        else
            YY = Ysyn; LL = double(Lsyn); 
        end
        LL(~LL)=-1;
    else
        YY = Y; LL = L; 
        if cfl, CC = C; end
        fl = false;
    end
    
end