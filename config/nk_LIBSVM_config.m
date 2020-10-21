function param = nk_LIBSVM_config(res, param, defaultsfl, cvfl, parentstr)

if ~exist('defaultsfl','var') || isempty(defaultsfl), defaultsfl=0; end
if ~exist('cvfl','var') || isempty(cvfl), cvfl = 0; end
if ~exist('param','var') || ~isfield(param,'LIBSVM'), param.LIBSVM = []; end

% DEFAULTS 
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% LIBSVM version
if isfield(param.LIBSVM,'LIBSVMver') , 
    LIBSVMver = param.LIBSVM.LIBSVMver;     
else
    LIBSVMver = 0;
end
switch LIBSVMver
    case 0
        libsvmstr = '3.12';
        quiet = '';
    case 1
        libsvmstr = '2.9.1';
        quiet = '';
    case 2
        libsvmstr = '2.89';
        quiet = ' -q 1';
    case 3
        libsvmstr = '2.89 PLUS (by Daniel Russo)';
        quiet = ' -q 1';
end
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Classifier / Regressor type
if isfield(param.LIBSVM,'classifier') && isnumeric(param.LIBSVM.classifier),
    classifier = param.LIBSVM.classifier;
else
    switch res.modeflag
        case 'classification'
            classifier = 0; 
        case 'regression'
            classifier = 3; 
    end
end
switch classifier
    case 0
        clstr = 'C-SVC (L1-regul.)'; 
    case 1
        switch LIBSVMver
            case {0,1,2}
                clstr = 'nu-SVC';
            case 3
                clstr = 'C-SVC (L2-regul.)';
        end
    case 2
        switch LIBSVMver
            case {0,1,2}
                clstr = 'one-class SVM';
            case 3
                clstr = 'nu-SVC';
        end
    case 3
        switch LIBSVMver
            case {0,1,2}
                clstr = 'epsilon-SVR';
            case 3
                clstr = 'one-class SVM';
        end
    case 4
        switch LIBSVMver
            case {0,1,2}
                clstr = 'nu-SVR';
            case 3
                clstr = 'epsilon-SVR';
        end
    case 5
        clstr = 'nu-SVR';
    case 6 
        clstr = 'SVDD (L1-regul.)';
    case 7
        clstr = 'SVDD (L2-regul.)';
end    
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Probability flag (Platt's method)
if isfield(param.LIBSVM,'Optimization') && isfield(param.LIBSVM.Optimization,'b') && isnumeric(param.LIBSVM.Optimization.b)
    probflag = param.LIBSVM.Optimization.b;
else
    probflag = 0;
end
switch probflag
    case 1
        probstr = 'yes';
    case 0
        probstr = 'no';
end
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Cachesize in MB
if isfield(param.LIBSVM,'Optimization') && isfield(param.LIBSVM.Optimization,'m') && isnumeric(param.LIBSVM.Optimization.m)
    cachesize = param.LIBSVM.Optimization.m;
else
    cachesize = 500;
end
cachestr = sprintf('%g MB',cachesize);
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Termination criterion
if isfield(param.LIBSVM,'Optimization') && isfield(param.LIBSVM.Optimization,'e') && isnumeric(param.LIBSVM.Optimization.e)
    termcrit = param.LIBSVM.Optimization.e;
else
    termcrit = 0.001;
end
termstr = sprintf('%g',termcrit);
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Shrinking heuristics
if isfield(param.LIBSVM,'Optimization') && isfield(param.LIBSVM.Optimization,'h') && isnumeric(param.LIBSVM.Optimization.h)
    shrinkheur = param.LIBSVM.Optimization.h;
else
    shrinkheur = 1;
end
switch shrinkheur
    case 1
        shrinkstr = 'yes';
    case 0
        shrinkstr = 'no';
end
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Weighting of hyperplane
if isfield(param.LIBSVM,'Weighting')
    weighting = param.LIBSVM.Weighting;
else
    weighting = 0;
end
switch weighting
    case 1
        weightstr = 'yes';
    case 0
        weightstr = 'no';
end

if isfield(param.LIBSVM,'W')
    caseweightstr = 'Weight vector defined';
else
    caseweightstr = 'Case weighting not activated';
end

% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Nu-Parameter
if isfield(param.LIBSVM,'Optimization') && isfield(param.LIBSVM.Optimization,'nu')
    nu = param.LIBSVM.Optimization.nu;
else
    nu = 0.5;
end
nustr = sprintf('%g',nu);
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Epsilon-Parameter
if isfield(param.LIBSVM,'Optimization') && isfield(param.LIBSVM.Optimization,'p')
    p = param.LIBSVM.Optimization.p;
else
    p = 0.8;
end
pstr = sprintf('%g',p);
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% INTERACTIVE LIBSVM MENU CONFIGURATION
if ~defaultsfl
    
    switch res.modeflag
        case 'classification'
            switch classifier 
                case 3
                    menuact = [ 'LIBSVM version [ ' libsvmstr ' ]|' ...
                                'Classifier type [ ' clstr ' ]|' ...
                                'Cache size [ ' cachestr ' ]|' ...
                                'Termination criterion [ ' termstr ' ]|' ...
                                'Shrinking heuristics [ ' shrinkstr ' ]|' ...
                                'Weighting of hyperplane in uneven group sizes [ ' weightstr ' ]|' ...
                                'Case-level weighting vector [ ' caseweightstr ' ]'];
                    menusel = [1,2,4:8];
                otherwise
                    menuact = [ 'LIBSVM version [ ' libsvmstr ' ]|' ...
                                'Classifier type [ ' clstr ' ]|' ...
                                'Probability estimation using Platt''s method [ ' probstr ' ]|' ...
                                'Cache size [ ' cachestr ' ]|' ...
                                'Termination criterion [ ' termstr ' ]|' ...
                                'Shrinking heuristics [ ' shrinkstr ' ]|' ...
                                'Weighting of hyperplane in uneven group sizes [ ' weightstr ' ]|' ...
                                'Case-level weighting vector [ ' caseweightstr ' ]'];
                    menusel = 1:8;
            end
        case 'regression'
             
             LIBSVMver = 0;
             menuact = ['Regressor type [ ' clstr ' ]|' ...
                        'Cache size [ ' cachestr ' ]|' ...
                        'Termination criterion [ ' termstr ' ]|' ...
                        'Shrinking heuristics [ ' shrinkstr ' ]|' ...
                        'Automatic case weighting according to target histogram [ ' weightstr ' ]|' ];
             menusel = [2,4:7];       
             
    end
    nk_PrintLogo
    mestr = 'LIBSVM: Parameter setup' ; navistr = [parentstr ' >>> ' mestr]; cprintf('*blue','\nYou are here: %s >>>',parentstr);
    act = nk_input( mestr,0,'mq', menuact, menusel);
    switch act
        case 1
            
           LIBSVMver = nk_input('LIBSVM version to use',0,'m',...
               ['LIBSVM 3.12 (supports data instance weighting)|' ...
               'LIBSVM 2.9.1|' ...
               'LIBSVM 2.89|' ...
               'LIBSVM 2.89 PLUS (only for classification)'], 0:3, LIBSVMver);
           switch LIBSVMver
               case {2,3}
                    quiet = ' -q 1';
               otherwise
                    quiet = '';
            end
                    
        case 2
            switch res.modeflag
                case 'classification'
                    switch LIBSVMver
                        case 3
                            menuclact = ['C-SVC (L1-regul.)|' ...
                                        'C-SVC (L2-regul.)|' ...
                                        'nu-SVC|' ...
                                        'one-class SVC|' ...
                                        'SVDD (L1-regul.)|' ...
                                        'SVDD (L2-regul.)|' ];
                            menuclsel = [0, 1, 2, 3, 6, 7];
                            cldef = find(menuclsel == classifier);
                        case {0,1,2}
                            menuclact = ['C-SVC (L1-regul.)|' ...
                                         'nu-SVC|' ...
                                         'one-class SVM'];
                            menuclsel = [0, 1, 2];
                            cldef = find(menuclsel == classifier);
                    end   
                case 'regression'
                    menuclact = ['epsilon-SVR (' pstr ')|' ...
                                'nu-SVR (' nustr ')'];
                    menuclsel = [3, 4];
                    cldef = find(menuclsel == classifier);     
            end
            classifier = nk_input('Classifier type to use',0,'m', menuclact, menuclsel, cldef);
            if ~cvfl && strcmp(res.modeflag,'regression')
                switch classifier
                    case 3
                        p = .1;
                    case 4
                        nu = .5;
                end
            end
        case 3
            if probflag, probflag = 0; else, probflag = 1; end
        case 4
            cachesize = nk_input('Cache size in MB',0,'e',cachesize);
        case 5
            termcrit = nk_input('Termination criterion',0,'e',termcrit);
        case 6
            if shrinkheur, shrinkheur = 0; else, shrinkheur = 1; end
        case 7
            if weighting , weighting = 0; else, weighting = 1; end
        case 8
            weightfl = nk_input('Enable case-level weighting',0,'y/n',[1,0]);
            if weightfl
                param.LIBSVM.W = nk_input('Vector with weights',0,'e',[],[numel(res.cases),1]);
            else
                if isfield(param.LIBSVM,'W')
                    param.LIBSVM = rmfield(param.LIBSVM,'W');
                end
            end
    end
else
    act = 0;
end

% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% TRANSFER PARAMS TO OUTPUT
param.LIBSVM.LIBSVMver = LIBSVMver;
param.LIBSVM.classifier = classifier;
param.LIBSVM.Optimization.b = probflag;
param.LIBSVM.Optimization.m = cachesize;
param.LIBSVM.Optimization.e = termcrit;
param.LIBSVM.Optimization.h = shrinkheur;
param.LIBSVM.Optimization.nu = nu;
param.LIBSVM.Optimization.p = p;
param.LIBSVM.Weighting = weighting;
param.LIBSVM.quiet = quiet;

if act, param = nk_LIBSVM_config(res, param, [], [], parentstr); end

