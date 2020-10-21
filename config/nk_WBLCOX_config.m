function [act, WBLCOX ] = nk_WBLCOX_config(NM, WBLCOX, defaultsfl)

if ~exist('defaultsfl','var') || isempty(defaultsfl), defaultsfl = false; end

if ~defaultsfl 
    if ~isfield(WBLCOX,'interval'),     WBLCOX.interval     = median(NM.time); end
    if ~isfield(WBLCOX,'MaxIter'),      WBLCOX.MaxIter      = 10000; end
    if ~isfield(WBLCOX,'MaxFunEvals'),  WBLCOX.MaxFunEvals  = 10000; end
    if ~isfield(WBLCOX,'TolX'),         WBLCOX.TolX         = 1e-12; end
    if ~isfield(WBLCOX,'TolFun'),       WBLCOX.TolFun       = 1e-12; end
    if ~isfield(WBLCOX,'Algorithm'),    WBLCOX.Algorithm    = 'sqp'; end
    if ~isfield(WBLCOX,'predict_time'), WBLCOX.predict_time = false; end
    PREDTIMEDEF = 'yes'; if ~WBLCOX.predict_time, PREDTIMEDEF = 'no'; end
    Algorithms = {'interior-point','trust-region-reflective','sqp','active-set'};
    AlgorithmDef = find(strcmp(Algorithms, WBLCOX.Algorithm));
    MnuStr = sprintf('Define time interval for risk estimation [ %g days ]|', WBLCOX.interval);                                 MnuAct = 1;
    MnuStr = sprintf('%sPredict time to event [ %s ]|', MnuStr, PREDTIMEDEF);                                                   MnuAct = [MnuAct 2];  
    MnuStr = sprintf('%sMaximum number of iterations allowed (positive scalar) [ %g ]|', MnuStr, WBLCOX.MaxIter);               MnuAct = [MnuAct 3];  
    MnuStr = sprintf('%sMaximum number of function evaluations allowed (positive scalar) [ %g ]|', MnuStr, WBLCOX.MaxFunEvals); MnuAct = [MnuAct 4];  
    MnuStr = sprintf('%sTermination tolerance on X (positive scalar)[ %g ]|', MnuStr, WBLCOX.TolX);                             MnuAct = [MnuAct 5];  
    MnuStr = sprintf('%sTermination tolerance on the function value (positive scalar) [ %g ]|', MnuStr, WBLCOX.TolFun);         MnuAct = [MnuAct 6];  
    MnuStr = sprintf('%sAlgorithm (menu-choice) [ %s ]', MnuStr, WBLCOX.Algorithm);                                             MnuAct = [MnuAct 7];
    
    nk_PrintLogo
    act = nk_input('Define Weibull-Cox Proportional Hazards Model parameters',0,'mq', MnuStr, MnuAct, 1);
    switch act
        case 1
            WBLCOX.interval     = nk_input('Define follow-up interval [days] for risk estimation',0,'e',WBLCOX.interval);
        case 2
            WBLCOX.predict_time = ~WBLCOX.predict_time;
        case 3
            WBLCOX.MaxIter = nk_input('Define number of algorithm iterations',0,'i',WBLCOX.MaxIter); 
        case 4
            WBLCOX.MaxFunEvals = nk_input('Define number of function evaluations',0,'i',WBLCOX.MaxFunEvals); 
        case 5
            WBLCOX.TolX = nk_input('Define tolerance of algorithm',0,'e',WBLCOX.TolX);
        case 6
            WBLCOX.TolFun = nk_input('Define tolerance of function evaluation',0,'e',WBLCOX.TolFun);
        case 7
            Algorithm = nk_input('Select algorithm',0,'mq',strjoin(Algorithms,'|'),1:numel(Algorithms),AlgorithmDef);
            WBLCOX.Algorithm = Algorithms{Algorithm};
    end
else
    act =0; WBLCOX = struct('interval',         median(NM.time), ...
                            'predict_time',     false, ...
                            'MaxIter',          1000, ...
                            'MaxFunEvals',      1000, ...
                            'TolX',             1e-05, ... 
                            'TolFun',           1e-05, ...
                            'Algorithm',        'sqp');
end


