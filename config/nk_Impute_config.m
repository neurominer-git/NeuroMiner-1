function [ IMPUTE, act ] = nk_Impute_config(NM, IMPUTE, varind, parentstr, defaultsfl)

IMPUTEDEF_K = 7; IMPUTEDEF_METH = 'euclidean' ; IMPUTEDEF_BLOCK = [];
if ~exist('varind','var') || isempty(varind), varind=1; end
IMPUTEDEF_HYBMETH1 = 'jaccard'; 
IMPUTEDEF_HYBMETH2 = 'euclidean';
IMPUTEDEF_HYBTHRESH = 3;

if ~exist('defaultsfl','var') || isempty(defaultsfl); defaultsfl = false; end

if ~defaultsfl
    
    if isfield(IMPUTE,'method'),    IMPUTEDEF_METH = IMPUTE.method;     end
    if isfield(IMPUTE,'blockind'),  IMPUTEDEF_BLOCK = IMPUTE.blockind; 	end
    if isfield(IMPUTE,'k'),         IMPUTEDEF_K = IMPUTE.k;             end
    if isfield(IMPUTE,'hybrid'),    
        IMPUTEDEF_HYBMETH1 = IMPUTE.hybrid.method1;
        IMPUTEDEF_HYBMETH2 = IMPUTE.hybrid.method2;
        IMPUTEDEF_HYBTHRESH = IMPUTE.hybrid.cutoff;
    end
    
    if isempty(IMPUTEDEF_BLOCK)
        IMPUTESTR_BLOCK = 'All features'; IMPUTEDEF_ALL = 1;
    else
        IMPUTESTR_BLOCK = sprintf('%g features',sum(IMPUTEDEF_BLOCK)); IMPUTEDEF_ALL = 2;
    end
   
    [ IMPUTESTR_METH, IMPUTEDEF_METH2 ] = return_methstr(IMPUTEDEF_METH);
    
    menustr = [ 'Define imputation method [ ' IMPUTESTR_METH ' ]|' ];
    menuact = 1;
    menustr = [ menustr 'Select features for imputation [ ' IMPUTESTR_BLOCK ' ]' ]; 
    menuact = [ menuact 2 ];
    
    switch IMPUTEDEF_METH2
         case {2,3,4,5,6,7,8} 
            IMPUTESTR_K = sprintf('%g nearest neighbors',IMPUTEDEF_K); 
            menustr = [menustr '|Define # of nearest neighbors [ ' IMPUTESTR_K ' ]'];
            menuact = [menuact 3];
            if IMPUTEDEF_METH2 == 8
                [ IMPUTESTR_HYBMETH1, IMPUTEDEF_HYBMETH_1 ] = return_methstr(IMPUTEDEF_HYBMETH1);
                [ IMPUTESTR_HYBMETH2, IMPUTEDEF_HYBMETH_2 ] = return_methstr(IMPUTEDEF_HYBMETH2);
                menustr = [menustr '|Define method for nominal features [ ' IMPUTESTR_HYBMETH1 ' ]']; menuact = [menuact 4];
                menustr = [menustr '|Define method for ordinal / continuous data [ ' IMPUTESTR_HYBMETH2 ' ]']; menuact = [menuact 5];
                menustr = sprintf('%s|Define maximum number of unique values for nominal feature imputation [ %g ]', menustr, IMPUTEDEF_HYBTHRESH); menuact = [menuact 6];
            end
    end
    
    nk_PrintLogo
    mestr = 'Imputation setup'; navistr = [parentstr ' >>> ' mestr]; cprintf('*blue','\nYou are here: %s >>> ',parentstr); 
    act = nk_input(mestr,0,'mq', menustr, menuact);
    
    switch act
        case 1
            IMPUTEDEF_METH = return_methdef(2, IMPUTEDEF_METH2);
        case 2
            allfeats = nk_input('Use all available features or define feature subspace',0,'m','All features|Subspace',[0,1], IMPUTEDEF_ALL);
            if allfeats
                if isfield(NM,'featnames') && ~isempty(NM.featnames{varind})
                    F = NM.featnames{varind};
                    fprintf('\n'); cprintf('*black','Feature list');
                    fprintf('\n'); cprintf('*black','============');
                    for i=1:numel(F)
                        fprintf('\n'); cprintf('*black',' [%4g]',i); fprintf('\t%s',F{i});
                    end   
                end
                blockind = nk_input('Define logical index of features to be imputed ',0,'e',[],size(NM.Y{varind},2));
                if ~islogical(blockind), blockind = logical(blockind); end
                IMPUTEDEF_BLOCK = false(1,size(NM.Y{varind},2));
                IMPUTEDEF_BLOCK(blockind) = true;
            else
                IMPUTEDEF_BLOCK = [];
            end
        case 3
            IMPUTEDEF_K = nk_input('# of nearest neighbors for computing the median',0,'w1',IMPUTEDEF_K);
        case 4
            IMPUTEDEF_HYBMETH1 = return_methdef(0, IMPUTEDEF_HYBMETH_1);
        case 5
            IMPUTEDEF_HYBMETH2 = return_methdef(1, IMPUTEDEF_HYBMETH_2);
        case 6
            IMPUTEDEF_HYBTHRESH = nk_input('Maximum number of unique values for nominal feature imputation',0,'w1',IMPUTEDEF_HYBTHRESH);
    end   
end

IMPUTE.blockind = IMPUTEDEF_BLOCK;
IMPUTE.k        = IMPUTEDEF_K; 
IMPUTE.method   = IMPUTEDEF_METH;
IMPUTE.hybrid.method1 = IMPUTEDEF_HYBMETH1;
IMPUTE.hybrid.method2 = IMPUTEDEF_HYBMETH2;
IMPUTE.hybrid.cutoff  = IMPUTEDEF_HYBTHRESH;
%==========================================================================

function [ IMPUTESTR_METH, IMPUTEDEF_METH2 ] = return_methstr(IMPUTEDEF_METH)

switch IMPUTEDEF_METH
    case 'singlemean'
        IMPUTESTR_METH = 'Replacement with SUBJECT-level mean across selected features'; IMPUTEDEF_METH2 = 0;
    case 'mean'
        IMPUTESTR_METH = 'GROUP-based mean imputation'; IMPUTEDEF_METH2 = 1;
    case 'cityblock'
        IMPUTESTR_METH = sprintf('kNN imputation (MANHATTAN)');  IMPUTEDEF_METH2 = 2;
    case 'euclidean'
        IMPUTESTR_METH = sprintf('kNN imputation (EUCLIDEAN)'); IMPUTEDEF_METH2 = 3;
    case 'seuclidean'
        IMPUTESTR_METH = sprintf('kNN imputation (Standardized EUCLIDEAN)'); IMPUTEDEF_METH2 = 4;
    case 'cosine'
        IMPUTESTR_METH = sprintf('kNN imputation (COSINE)'); IMPUTEDEF_METH2 = 5;
    case 'hamming'
        IMPUTESTR_METH = sprintf('kNN imputation (HAMMING: nominal only)'); IMPUTEDEF_METH2 = 6;
    case 'jaccard'
        IMPUTESTR_METH = sprintf('kNN imputation (JACCARD: nominal positive only)'); IMPUTEDEF_METH2 = 7;
    case 'hybrid'
        IMPUTESTR_METH = sprintf('kNN imputation (HYBRID: Adaptive to nominal vs. ordinal/continuous features'); IMPUTEDEF_METH2 = 8;
end

function IMPUTEDEF_METH = return_methdef(TYPE, IMPUTEDEF_METH2)

switch TYPE 
    case 0 % Nominal data
        if exist('pdist2','file')
            IMPUTEDEF_METHMENU = ['HAMMING distance-based nearest-neighbor search', ...
                                  '|JACCARD distance-based nearest-neighbor search'];
            IMPUTEDEF_METHMENUSEL = {'hamming', 'jaccard'};
       end
    
    case 1 % Coninuous data
        IMPUTEDEF_METHMENU = ['MANHATTAN distance-based nearest-neighbor search (Scale / Standardize first!)|', ...
                 'EUCLIDEAN distance-based nearest-neighbor search (Scale / Standardize first!)'];
        if ~isempty(which('pdist2'))
            IMPUTEDEF_METHMENUSEL = {'cityblock', 'euclidean'};
            IMPUTEDEF_METHMENU = [IMPUTEDEF_METHMENU ...
                                  '|SEUCLIDEAN distance-based nearest-neighbor search', ...
                                  '|COSINE similarity-based nearest-neighbor search'];
            IMPUTEDEF_METHMENUSEL = [IMPUTEDEF_METHMENUSEL  {'seuclidean', 'cosine'}];
        else
            error('Statistics toolbox not available!')
        end
    case 2 % All
        IMPUTEDEF_METHMENU = ['Median of non-NaN values in given case|', ...
                 'Mean of non-NaN values in given feature|', ...
                 'MANHATTAN distance-based nearest-neighbor search (Scale / Standardize first!)|', ...
                 'EUCLIDEAN distance-based nearest-neighbor search (Scale / Standardize first!)'];
        IMPUTEDEF_METHMENUSEL = {'singlemean', 'mean', 'cityblock', 'euclidean'};
        if  ~isempty(which('pdist2'))
            IMPUTEDEF_METHMENU = [IMPUTEDEF_METHMENU ...
                                  '|SEUCLIDEAN distance-based nearest-neighbor search', ...
                                  '|COSINE similarity-based nearest-neighbor search', ...
                                  '|HAMMING distance-based nearest-neighbor search', ...
                                  '|JACCARD distance-based nearest-neighbor search', ...
                                  '|Nearest-neighbor imputation using hybrid method'];
            IMPUTEDEF_METHMENUSEL = [IMPUTEDEF_METHMENUSEL  {'seuclidean', 'cosine', 'hamming', 'jaccard', 'hybrid'}];
        else
             error('Statistics toolbox not available!')
        end
end
IMPUTEDEF_METH = char(nk_input('Replacement of NaN by',0,'m', IMPUTEDEF_METHMENU,  IMPUTEDEF_METHMENUSEL , IMPUTEDEF_METH2));