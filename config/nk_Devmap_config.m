function [DEVMAP, PX] = nk_Devmap_config(DEVMAP, PX, NM, defaultsfl, parentstr)

if ~exist('defaultsfl','var') || isempty(defaultsfl),  defaultsfl = 0; end;
algostr     = 'pls';
covdesc     = 'NM target label';
covmat      = NM.label;
loo         = false;
glabel      = [];
ncomp       = 1;
act         = 0;

if ~defaultsfl

    if isfield(DEVMAP,'algostr'), algostr = DEVMAP.algostr; end
    if isfield(DEVMAP,'covmat'),  covmat = DEVMAP.covmat; end
    if isfield(DEVMAP,'covdesc'), covdesc = DEVMAP.covdesc; end
    if isfield(DEVMAP,'glabel'),  glabel = DEVMAP.glabel; end
    if isfield(DEVMAP,'ncomp'),   ncomp = DEVMAP.ncomp; else, ncomp = size(covmat,2); end
    if isfield(DEVMAP,'PX'),      PX = DEVMAP.PX; end
    if isfield(DEVMAP,'cu'),      cu = DEVMAP.cu; end
    if isfield(DEVMAP,'cv'),      cv = DEVMAP.cv; end
    nk_PrintLogo
    
    if ~isempty(glabel)
        grpstr = sprintf('yes, compute model in %g cases', sum(glabel));
    else
        grpstr = sprintf('no');
    end
    
    if size(covmat,2)>1
        covsizestr = sprintf('Matrix ''%s'' (%g x %g)', covdesc, size(covmat,1), size(covmat,2));
        componentstr = []; mnuact = 1:3;
    else
        covsizestr = sprintf('Vector ''%s'' (%g x 1)', covdesc, size(covmat,1));
        componentstr = []; mnuact = 1:3;
    end
    if strcmp(algostr,'spls')
        if ~exist('cu','var'), yyy = 'undefined'; cudef = 1; else, yyy = nk_ConcatParamstr(cu); cudef = cu; end 
        if ~exist('cv','var'), xxx = 'undefined'; cvdef = 1; else, xxx = nk_ConcatParamstr(cv); cvdef = cv; end
        sparsepls = sprintf('|Define U sparsity parameter range (1 to sqrt(input dimensionality)) [ %s ]', yyy);
        mnuact = [mnuact 4];
        if size(covmat,2)>1
            sparsepls = sprintf('%s|Define V sparsity parameter range (1 to sqrt(input dimensionality)) [ %s ]', sparsepls, xxx);
            mnuact = [mnuact 5];
        else
            PX = nk_AddParam(1, 'SPLS-cv',1, PX);
        end
    else
        sparsepls=[];
    end
    mestr = 'Deviation-based feature weighting'; navistr = [parentstr ' >>> ' mestr]; cprintf('*blue','\nYou are here: %s >>> ',parentstr); 
    
    act = nk_input(mestr,0,'mq', ...
        [sprintf('Choose algorithm and specify its parameters [ %s ]|', algostr), ...
         sprintf('Compute normative model in a subgroup [ %s ]|', grpstr), ...
         sprintf('Define covariate matrix for PLS-based deviation ranking [ %s ]', covsizestr), ...
         componentstr ...
         sparsepls], mnuact);
    
    switch act
        case 1
            algostr = char(nk_input('Select normative algorithm',0,'mq','Partial Least Squares|Sparse Partial Least Squares',{'pls','spls'},1));
            if strcmp(algostr,'pls'), PX = nk_AddParam([], [], [], PX,'reset'); end
        case 2
            flg = nk_input('Weight feature using only one specific subgroup?',0,'yes|no',[1,0],0);
            if flg, glabel = nk_input('Define logical vector to identify subgroup',0,'e',[],[numel(NM.label),1]); else, glabel=[]; end
        case 3
            covmat = nk_input(sprintf('Define covariate matrix for %s',algostr),0,'e',[],[numel(NM.label),Inf]); 
            covdesc = nk_input('Give a short description for the covariate matrix/vector',0,'s');
            ncomp = size(covmat,2);
        case 4
            DEVMAP.cu = nk_input('Define U sparsity parameter(s)',0,'e',cudef);
            PX = nk_AddParam(DEVMAP.cu, 'SPLS-cu',1, PX);
        case 5
            DEVMAP.cv = nk_input('Define V sparsity parameter(s)',0,'e',cvdef);
            PX = nk_AddParam(DEVMAP.cv, 'SPLS-cv',1, PX);
    end
end
DEVMAP.algostr      = algostr;
DEVMAP.covmat       = covmat;
DEVMAP.covdesc      = covdesc;
DEVMAP.glabel       = glabel;
DEVMAP.ncomp        = ncomp;
DEVMAP.loo          = loo;
DEVMAP.PX           = PX;
if act, [DEVMAP, PX] = nk_Devmap_config(DEVMAP, PX, NM, defaultsfl, parentstr); end
