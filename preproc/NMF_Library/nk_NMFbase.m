function [ w_sol, infos_sol, costs ] = nk_NMFbase( V, rank, algo, options)

[F, N] = size(V);
Vo = V;
dim = N * F;

%% set initial data
x_init.W = rand(F, rank); 
x_init.H = rand(rank, N);
x_init.R = rand(F, N);

%% set options
if ~exist('options', 'var') || isempty(options)
	max_epoch = 100;	
    options.max_epoch = max_epoch;
    options.x_init = x_init;
    options.verbose = true;  
    options.max_epoch = 1000;
	options.calc_sol = false;
    if strcmp(algo,'nenmf')
        options.type = 'l1r';
        options.stop_rule = 1;
    end
end

%% calculate optimal solution
if options.calc_sol
    fprintf('Calculating f_opt by HALS ...\n');
    options.alg = 'hals';
    w_sol = nmf_als(V, rank, options);
    f_opt = nmf_cost(V, w_sol.W, w_sol.H, zeros(F, N));    
    fprintf('Done.. f_opt: %.16e\n', f_opt);
else
    f_opt = -Inf;
end

options.f_opt = f_opt;   
options.alg = algo;
algostr = [ 'nmf_' algo ];

switch algo
	case 'mu_mod'
		algostr = 'nmf_mu';
	case 'mu_acc'
		options.alpha = 2;
		options.delta = 0.1;
	case {'hals','hals_acc'}
		algostr = 'nmf_als';
	case 'direct_pgd'
		algostr = 'nmf_pgd';
	case  'anls_bpp'
		algostr = 'nmf_anls';
    case 'nenmf'
        algostr = algo;
end
[w_sol, infos_sol] 	= feval(algostr, V, rank, options);
costs = nmf_cost(Vo, w_sol.W, w_sol.H, zeros(F, N)) * 2 / dim;
