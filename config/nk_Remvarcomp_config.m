function [REMVARCOMP, PX] = nk_Remvarcomp_config(NM, varind, REMVARCOMP, PX, parentstr)

[m, ~] = size(NM.Y{varind});
corrmeth    = 1;
corrthresh  = 0.5;
varop = 'gt';
recon = 2;
varops = {'gt','ge','lt','le'};
recons = {'yes', 'no'};
corrmeths = {'Pearson','Spearman','ANOVA'};
dims = 0.9;
dimmode = 3;
SUBGROUP.flag = 1;
subgroupstr = {'Entire training set', ...
                'User-defined subset of CV1 training partition', ...
                'Random subset of CV1 training partition', ...
                'Entire calibration set', ...
                'User-defined subset of calibration data',...
                'Random subset of caobration data'};
            
%% Define defaults
if ~isfield(REMVARCOMP,'corrmeth'),     REMVARCOMP.corrmeth = corrmeth; end
if ~isfield(REMVARCOMP,'corrthresh'),   REMVARCOMP.corrthresh = corrthresh; end
if ~isfield(REMVARCOMP,'varop'),        REMVARCOMP.varop = varop; end
if ~isfield(REMVARCOMP,'recon'),        REMVARCOMP.recon = recon; end
if ~isfield(REMVARCOMP,'dims'),         REMVARCOMP.dims = dims; end
if ~isfield(REMVARCOMP,'dimmode'),      REMVARCOMP.dimmode = dimmode; end
if ~isfield(REMVARCOMP,'SUBGROUP'),     REMVARCOMP.SUBGROUP = SUBGROUP; end
if ~exist('PX','var'),                  PX = []; end

%% Define menu
if isfield(REMVARCOMP,'G')
    [p,q] = size(REMVARCOMP.G);
    if q>1
        RANKVARCOMP_G_str = sprintf('%g-by-%g matrix specified',p,q);
    else
        RANKVARCOMP_G_str = 'vector specified';
    end
else
    RANKVARCOMP_G_str = 'No covariate vector or matrix specified!';
end

REMVARCOMP_corrmeth_str = corrmeths{REMVARCOMP.corrmeth};
REMVARCOMP_corrthresh_str = nk_ConcatParamstr(REMVARCOMP.corrthresh);
varopdef = find(strcmp(varops, REMVARCOMP.varop));
REMVARCOMP_recon_str = recons{REMVARCOMP.recon};
switch REMVARCOMP.dimmode
    case 1
        REMVARCOMP_dims_str = ['Fixed: ' nk_ConcatParamstr(REMVARCOMP.dims) ' components '];
    case 3
        REMVARCOMP_dims_str = ['Percentage: ' nk_ConcatParamstr(REMVARCOMP.dims) ' of total variance '];
end
REMVARCOMP_subgroup_str = subgroupstr{REMVARCOMP.SUBGROUP.flag};

menustr = [ 'Define target vector / matrix for identification of variance components [ ' RANKVARCOMP_G_str ' ]|', ...
            'Define correlation method for identification of variance components [ ' REMVARCOMP_corrmeth_str ' ]|', ...
            'Define correlation cutoff for variance extraction [ ' REMVARCOMP_corrthresh_str  ' ]|', ...
            'Define variance extraction operator [ ' REMVARCOMP.varop  ' ]|', ...
            'Back-project adjusted matrices to input space [ ' REMVARCOMP_recon_str ' ]|'...
            'Define % retained variance in PCA projection [ ' REMVARCOMP_dims_str ']|' ...
            'Define training data for PCA model generation [ ' REMVARCOMP_subgroup_str ' ]'];
menuact = 1:7;

%% Print menu
nk_PrintLogo
mestr = 'Extract variance components from matrix'; navistr = [parentstr ' >>> ' mestr]; cprintf('*blue','\nYou are here: %s >>> ',parentstr); 
act = nk_input(mestr,0,'mq', menustr , menuact);

% Get user data
switch act
    case 1
        REMVARCOMP.G = nk_input('Define target vector (matrix) for identification of variance components',0,'e',[],[m Inf]);
    case 2
        REMVARCOMP.corrmeth = nk_input('Define method for identifying variance components',0,'m','Pearson|Spearman|ANOVA', [1,2,3], REMVARCOMP.corrmeth);
    case 3
        REMVARCOMP.corrthresh = nk_input('Define single absolute correlation cutoff or multiple cutoffs for variance extraction',0,'e',REMVARCOMP.corrthresh);
    case 4
        varopnum = nk_input('Define variance extraction operator',0,'>|>=|<|<=',{'gt','ge','lt','le'}, varopdef); REMVARCOMP.varop = char(varopnum);
    case 5
        if REMVARCOMP.recon == 1, REMVARCOMP.recon = 2; elseif REMVARCOMP.recon == 2, REMVARCOMP.recon = 1; end
    case 6
        switch REMVARCOMP.dimmode
            case 1
                dimmodedef = 1;
            case 3
                dimmodedef = 2;
        end
        REMVARCOMP.dimmode = nk_input('Retain fixed number or percentage of variance components',0,'m','Fixed|Percentage',[1,3],dimmodedef);
        switch REMVARCOMP.dimmode
            case 1
                REMVARCOMP.dims = nk_input('Fixed number of variance components to be returned by PCA (1-num(dims))',0,'e',REMVARCOMP.dims); 
            case 3
                REMVARCOMP.dims = nk_input('Percentage of variance components to be returned by PCA (0-1)',0,'e',REMVARCOMP.dims); 
        end
    case 7
        REMVARCOMP.SUBGROUP = nk_SubgroupOperation_config( NM, REMVARCOMP.SUBGROUP );
end

% Register correlation threshold to parameter space if numel(thresholds) > 1 
PX = nk_AddParam(REMVARCOMP.corrthresh, 'CorrThresh', 1, PX, 'replace');
PX = nk_AddParam(REMVARCOMP.dims, 'Dims', 1, PX); 
if act, [REMVARCOMP, PX] = nk_Remvarcomp_config(NM, varind, REMVARCOMP, PX, parentstr); end