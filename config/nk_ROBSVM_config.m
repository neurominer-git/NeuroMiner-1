function [ param, act ] = nk_ROBSVM_config(prog, param, modefl, defaultsfl)


if strcmp(modefl,'classification'), modefl = 'binaryclass'; end
if ~exist('defaultsfl','var') || isempty(defaultsfl), defaultsfl = true; end

probflag = 0;
cachesize = 500;
shrinkheur = 1;
termcrit = 0.001;
weighting = 0;
act = 0;
if isfield(param,'Optimization')
    probflag = param.Optimization.b;
    cachesize = param.Optimization.m; 
    shrinkheur = param.Optimization.h; 
    termcrit = param.Optimization.e; 
end
if isfield(param,'Weighting'); weighting = param.Weighting; end

param.algo = prog;
param.learner.framework = modefl;
param.learner.options = nk_matLearn_getopts_config([], 'get_learner_params', prog, param, modefl);
param = nk_matLearn_DefAlgoParams_config(param, param.learner.options, false);

if ~defaultsfl
    
    switch probflag
        case 1
            probstr = 'yes';
        case 0
            probstr = 'no';
    end
    cachestr = sprintf('%g MB',cachesize);
    termstr = sprintf('%g',termcrit);
    switch shrinkheur
        case 1
            shrinkstr = 'yes';
        case 0
            shrinkstr = 'no';
    end
    switch weighting
        case 1
            weightstr = 'yes';
        case 0
            weightstr = 'no';
    end
    nk_PrintLogo
    
    act = nk_input('Define ROBSVM settings',0,'mq',[ 'Define hyperparameters|' ...
                            'Define probabilistic output [' probstr ']|' ...
                            'Cache size [ ' cachestr ' ]|' ...
                            'Termination criterion [ ' termstr ' ]|' ...
                            'Shrinking heuristics [ ' shrinkstr ' ]|' ...
                            'Weighting of hyperplane for unbalanced problems [' weightstr ']'],1:6);

    switch act
        case 1
            param = nk_matLearn_DefAlgoParams_config(param, param.learner.options, true);
        case 2
            param.Optimization.b =  nk_input('Use probability estimates (this leads to slower computations)',0,'yes|no',[1,0],find([1 0] == probflag));
        case 3
            param.Optimization.m = nk_input('Cache size in MB',0,'e',cachesize);
        case 4
            param.Optimization.e = nk_input('Termination criterion',0,'e',termcrit);
        case 5
            param.Optimization.h = nk_input('Use shrinking heuristics',0,'yes|no',[1,0],find([1 0] == shrinkheur));
        case 6
            param.Weighting = nk_input('Weight hyperplane',0,'yes|no',[1,0],find([1 0] == weighting));
    end
        
else
    for ind = 1:numel(param.learner.options.name)
        param.Params(ind).range = param.learner.options.def{ind};
        param.Params(ind).name = param.learner.options.name{ind}; 
    end
    param.Optimization = struct('b',probflag, 'm',cachesize,'e', termcrit, 'h', shrinkheur);
    param.Weighting = weighting;
     
end


