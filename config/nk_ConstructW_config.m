function [Param, PX] = nk_ConstructW_config(Param, defaultsfl, parentstr)

if ~exist('defaultsfl','var') || isempty(defaultsfl), defaultsfl=0; end

Metric = 'Cosine';
NeighborMode = 'KNN';
WeightMode = 'Binary';
bNormalized = 0;
bLDA = 0;
K = 5;
T = 1;
PX = [];

if ~defaultsfl
    
    if isempty(Param)
        Param.Metric        = Metric;
        Param.NeighborMode  = NeighborMode;
        Param.WeightMode    = WeightMode;
        Param.bNormalized   = bNormalized;
        Param.bLDA          = bLDA;
        Param.k             = K;
        Param.t             = T;
    else
        Metric              = Param.Metric;
        NeighborMode        = Param.NeighborMode;
        WeightMode          = Param.WeightMode;
        bNormalized         = Param.bNormalized;
        bLDA                = Param.bLDA;
        K                   = Param.k;
        T                   = Param.t;
    end
    
    switch Param.Metric
        case 'Euclidean'
            MetricDef = 1; 
        case 'Cosine'
            MetricDef = 2;
    end
    switch Param.NeighborMode
        case 'KNN'
            NeighborModeDef = 1;
        case 'Supervised'
            NeighborModeDef = 2;
    end
    switch Param.WeightMode
        case 'Binary'
            WeightModeDef = 1;
        case {'HeatKernel','Cosine'}
            WeightModeDef  = 2;
    end
    switch Param.bNormalized
        case 1
            bNormalizedDef = 1; bNormalizedstr = 'yes';
        case 0
            bNormalizedDef = 2; bNormalizedstr = 'no';
    end
    switch Param.bLDA
        case 1
            bLDADef = 1; bLDAstr = 'yes';
        case 0
            bLDADef = 2; bLDAstr = 'no';
    end
    KDef = Param.k; TDef = T;
    
    menustr = sprintf('How shall we assess closeness of two datapoints ? [ %s ]', Param.Metric); menuact = 1;
    menustr = sprintf('%s|Define the method to be used for constructing the graph ? [ %s ]', menustr, Param.NeighborMode); menuact = [ menuact 2 ];
    menustr = sprintf('%s|How are weights assigned to the edges in the graph ? [ %s ]', menustr, Param.WeightMode); menuact = [ menuact 3 ];
    switch Param.Metric
        case 'Euclidean'
            weightmode_mnu  = 'Binary|HeatKernel';
            weightmode_sel  = {'Binary','HeatKernel'};
            if strcmp(Param.WeightMode,'HeatKernel')
                menustr         = sprintf('%s|Heat kernel parameter t [ %s ]', menustr, nk_ConcatParamstr(Param.t)); menuact = [ menuact 7];
                PX = nk_AddParam(Param.t, 'LPP-T', 1, PX);
            end
        case 'Cosine'
            weightmode_mnu  = 'Binary|Cosine'; 
            weightmode_sel  = {'Binary','Cosine'}; 
            menustr         = sprintf('%s|Are the features already normalized ? (Cosine metric only) [ %s ]', menustr, bNormalizedstr); menuact = [ menuact 4];
    end
    if strcmp(Param.NeighborMode,'Supervised')
        menustr = sprintf('%s|Use LDA to construct the graph ? (Supervised graph construction only) [ %s ]', menustr, bLDAstr); menuact = [ menuact 5];
    end
    menustr = sprintf('%s|How many nearest neighbors should be taken into account ? [ %s ]', menustr, nk_ConcatParamstr(Param.k)); menuact = [ menuact 6];
    PX = nk_AddParam(Param.k, 'LPP-K', 1, PX);
    
    nk_PrintLogo
    mestr = 'Graph construction (Affinity matrix) parameters'; navistr = [parentstr ' >>> ' mestr]; cprintf('*blue','\nYou are here: %s >>> ',parentstr); 
    act = nk_input(mestr,0,'mq', menustr, menuact);
    
    switch act
        case 1
            Metric        = char(nk_input('Affinity Matrix Setup: Metric',0,'m','Euclidean|Cosine',{'Euclidean','Cosine'}, MetricDef));
        case 2
            NeighborMode  = char(nk_input('Affinity Matrix Setup: NeighborMode',0,'m','KNN|Supervised',{'KNN','Supervised'},NeighborModeDef));
        case 3
            WeightMode    = char(nk_input('Affinity Matrix Setup: Weight Mode',0, weightmode_mnu, weightmode_sel, WeightModeDef ));
        case 4
            bNormalized   = nk_input('Affinity Matrix Setup: Normalized (Cosine)',0,'yes|no',[1,0], bNormalizedDef);
        case 5
            bLDA = nk_input('Affinity Matrix Setup: bLDA (Supervised Neighbormode)',0,'yes|no',[1,0], bLDADef);
        case 6
            K = nk_input('Affinity Matrix Setup: k (nearest neighbors)',0,'e', KDef); 
            PX = nk_AddParam(K, 'LPP-K', 1, PX, 'replace');
        case 7
            T = nk_input('Affinity Matrix Setup: t parameter range (Heatm Kernel param)',0,'e', TDef); 
            PX = nk_AddParam(T, 'LPP-T', 1, PX, 'replace');
    end
else
    act = 0;
end

Param.Metric        = Metric;
Param.NeighborMode  = NeighborMode;
Param.WeightMode    = WeightMode;
Param.bNormalized   = bNormalized;
Param.bLDA          = bLDA;
Param.k             = K;
Param.t             = T;

if act, [Param, PX] = nk_ConstructW_config(Param, [], parentstr); end