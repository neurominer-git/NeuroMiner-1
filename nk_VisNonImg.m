function nk_VisNonImg(inp, strout, id, GridAct, batchflag)

global SVM SAV DR RFE MULTI SCALE DISCRET SYMBOL COVAR FEATSEL CLUST MODEFL VIS

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%% INITIALIZATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
Y               = inp.Y;            % Input data matrix
nclass          = inp.nclass;       % # of binary comparisons
cv              = inp.cv;           % cross-validation structure
dim             = inp.dimension;    % current dimensionality
analysis        = inp.analysis;     % GDanalysis structure to be used
covstr          = inp.covstr;       % Covariate description

switch inp.analmode
    case 0
        labels          = inp.labels;       % Labels / Targets
        covars          = inp.covars;       % Covariates
        preprocmat      = inp.preprocmat;   % paths to preprocessed data files
        gdmat           = inp.gdmat;        % paths to precomputed GD structures
        featmat         = inp.featmat;      % paths to preprocessed feature selection files
        ovrwrt          = inp.ovrwrt;       % overwrite existing data files
    case 1
        vismat          = inp.vismat;       % Visualization datamat
end

% Check whether RVM or SVM will be used?
if strcmp(SVM.prog,'MikSVM'),
    if strcmp(MODEFL,'classification')
        algostr = 'RVM';
    else
        algostr = 'RVR';
    end
else
    if strcmp(MODEFL,'classification')
        algostr = 'SVM';
    else
        algostr = 'SVR';
    end
end

% Setup CV2 container variables:
[ix, jx] = size(cv.TrainInd);

% Get CV2gridflags if not specified
if ~exist('GridAct','var') || isempty(GridAct), ...
        GridAct=nk_CVGridSelector(ix,jx); end

if ~exist('batchflag','var') || isempty(batchflag), batchflag = false; end

ol = 0; binmode = FEATSEL.binmode;
VCV2SUM = zeros(size(Y,2),nclass);
VCV2SQ = zeros(size(Y,2),nclass);
GCV2SUM = nan(size(Y,2),ix*jx,nclass);
ind0 = false(1,ix*jx);
[iy, jy] = size(cv.cvin{1,1}.TrainInd);
ll = 1; numCV1partsG = zeros(nclass,1);

for f=1:ix % Loop through CV2 permutations

    for d=1:jx % Loop through CV2 folds
        
        fprintf('\n--------------------------------------------------------------------------')
        
    end
end