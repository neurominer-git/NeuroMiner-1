function param = nk_SA_config(param, defaultsfl)

SA.c                = 0.01; % Sparsity constant
SA.T                = 1;    % Starting temperature
SA.T_stop           = 0.005; % Stopping temperature 
SA.alpha            = 0.9;  
SA.itt_max          = 2000; % Max number of iterations
SA.Rep_T_max        = 100;  % max number of Rep in one temperature T
SA.Rep_accept_max   = 30;   % if solution accepted this many, then update T
SA.kc               = 10;   % k-constant. The smaller kc --> less solution accepted @#$% user-defined

if ~exist('defaultsfl','var') || isempty(defaultsfl), defaultsfl = false; end
if ~defaultsfl
    
    if isfield(param,'SA') 
        if isfield(param.SA,'c') && ~isempty( param.SA.c ), SA.c = param.SA.c; end
        if isfield(param.SA,'T') && ~isempty( param.SA.T ), SA.T = param.SA.T; end
        if isfield(param.SA,'T_stop') && ~isempty( param.SA.T_stop ), SA.T_stop = param.SA.T_stop; end
        if isfield(param.SA,'alpha') && ~isempty( param.SA.alpha ), SA.alpha = param.SA.alpha; end
        if isfield(param.SA,'itt_max') && ~isempty( param.SA.itt_max ), SA.itt_max = param.SA.itt_max; end
        if isfield(param.SA,'Rep_T_max') && ~isempty( param.SA.Rep_T_max ), SA.Rep_T_max = param.SA.Rep_T_max; end
        if isfield(param.SA,'Rep_accept_max') && ~isempty( param.SA.Rep_accept_max ), SA.Rep_accept_max = param.SA.Rep_accept_max; end
        if isfield(param.SA,'kc') && ~isempty( param.SA.kc ), SA.kc = param.SA.kc; end
    end
    nk_PrintLogo
    act = nk_input('Define Simulated annealing parameters',0,'mq', ...
                    [sprintf('Sparsity constant [ %g ]|', SA.c) ...
                     sprintf('Starting temperature [ %g ]|', SA.T) ...
                     sprintf('Stopping temperature [ %g ]|', SA.T_stop) ...
                     sprintf('Alpha [ %g ]|', SA.alpha) ...
                     sprintf('Maximum # of iterations [ %g ]|', SA.itt_max) ...
                     sprintf('Maximum # of repetitions in Temperature T [ %g ]|', SA.Rep_T_max) ...
                     sprintf('Minimum # of repetitions until solution is accepted in T [ %g ]|', SA.Rep_accept_max) ...
                     sprintf('k constant (smaller k ==> less solutions accepted) [ %g ]', SA.kc) ], 1:8);
    switch act
        case 1
            SA.c = nk_input('Define SA Sparsity Constant',0,'e',SA.c);
        case 2
            SA.T = nk_input('Define SA Starting Temperature',0,'e',SA.T);
        case 3
            SA.T_stop = nk_input('Define SA Stopping Temperature',0,'e',SA.T_stop);
        case 4
            SA.alpha = nk_input('Define SA Alpha',0,'e',SA.alpha);
        case 5
            SA.itt_max = nk_input('Define SA maximum # of iterations',0,'e',SA.itt_max);
        case 6
            SA.Rep_T_max = nk_input('Define maximum # of repetitions in Temperature T',0,'e',SA.Rep_T_max);
        case 7
            SA.Rep_accept_max = nk_input('Define minimum # of repetitions until solution is accepted in T',0,'e',SA.Rep_accept_max);
        case 8
            SA.kc = nk_input('Define SA k constant',0,'e',SA.kc);
    end
else
    act = 0;
end
param.SA = SA; 
if act, param = nk_SA_config(param, defaultsfl); end