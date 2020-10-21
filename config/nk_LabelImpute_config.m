function [ LABELIMPUTE, act ] = nk_LabelImpute_config(LABELIMPUTE, parentstr, defaultsfl)

IMPUTEDEF_K = 7; IMPUTEDEF_METH = 'none' ; 

if ~exist('defaultsfl','var') || isempty(defaultsfl); defaultsfl = false; end

if ~defaultsfl
    
    if isfield(LABELIMPUTE,'method'),    IMPUTEDEF_METH = LABELIMPUTE.method;     end
    if isfield(LABELIMPUTE,'k'),         IMPUTEDEF_K    = LABELIMPUTE.k;          end
        
    switch IMPUTEDEF_METH
        case 'none'
            IMPUTESTR_METH = sprintf('No label propagation');                           IMPUTEDEF_METH2 = 1;
        case 'ml'
            IMPUTESTR_METH = sprintf('Label propagation using the supervised learner'); IMPUTEDEF_METH2 = 2;
        case {'manhattan'}
            IMPUTESTR_METH = sprintf('MANHATTAN-based KNN replacement');                IMPUTEDEF_METH2 = 3;
        case {'euclidean'}
            IMPUTESTR_METH = sprintf('EUCLIDEAN-based KNN replacement');                IMPUTEDEF_METH2 = 4;
        case 'seuclidean'
            IMPUTESTR_METH = sprintf('Standardized EUCLIDEAN-based KNN replacement');   IMPUTEDEF_METH2 = 5;
        case 'cosine'
            IMPUTESTR_METH = sprintf('COSINE-based KNN replacement');                   IMPUTEDEF_METH2 = 6;
        case 'hamming'
            IMPUTESTR_METH = sprintf('HAMMING-based KNN replacement');                  IMPUTEDEF_METH2 = 7;
        case 'jaccard'
            IMPUTESTR_METH = sprintf('JACCARD-based KNN replacement');                  IMPUTEDEF_METH2 = 8;
    end
    
    menustr = [ 'Define imputation method [ ' IMPUTESTR_METH ' ]|' ];
    menuact = 1;
    
    if IMPUTEDEF_METH2>1
        IMPUTESTR_K = sprintf('%g nearest neighbors',IMPUTEDEF_K); 
        menustr = [menustr '|Define # of nearest neighbors [ ' IMPUTESTR_K ' ]'];
        menuact = [menuact 2];
    end
    
    nk_PrintLogo
    mestr = 'Imputation setup'; navistr = [parentstr ' >>> ' mestr]; cprintf('*blue','\nYou are here: %s >>> ',parentstr); 
    act = nk_input(mestr,0,'mq', menustr, menuact);
    
    switch act
        case 1
            IMPUTEDEF_METHMENU = [    'No label propagation', ...
                                      '|Label propagation using the selected supervised machine learning algorithm', ... 
                                      '|MANHATTAN distance-based nearest-neighbor search (Scale / Standardize first!)', ...
                                      '|EUCLIDEAN distance-based nearest-neighbor search (Scale / Standardize first!)', ...
                                      '|SEUCLIDEAN distance-based nearest-neighbor search', ...
                                      '|COSINE similarity-based nearest-neighbor search', ...
                                      '|HAMMING distance-based nearest-neighbor search', ...
                                      '|JACCARD distance-based nearest-neighbor search' ];
                                  
            IMPUTEDEF_METHMENUSEL = {'none','ml','manhattan', 'euclidean', 'seuclidean', 'cosine', 'hamming', 'jaccard'};
            
            IMPUTEDEF_METH = char(nk_input('Replacement of NaN by',0,'m', IMPUTEDEF_METHMENU,  IMPUTEDEF_METHMENUSEL , IMPUTEDEF_METH2));
        case 2
            IMPUTEDEF_K = nk_input('# of nearest neighbors for computing the mean',0,'w1',IMPUTEDEF_K);
    end   
end
%==========================================================================
LABELIMPUTE.k        = IMPUTEDEF_K; 
LABELIMPUTE.method   = IMPUTEDEF_METH;