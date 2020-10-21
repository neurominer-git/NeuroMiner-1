% ==========================================================================
function opt = nk_matLearn_getopts_config(opt, act, algo, param, framework)

if ~exist('param','var'), param = []; end

switch act
    
    case {'get_learners','get_sublearners'}
        
        opt.desc    = sprintf('Select %s algorithm from the matLearn library',framework);
        opt.name    = 'algo';
        opt.def     = 1;
        
        switch act
            case 'get_learners'
                switch framework
                    case 'regression'
                            opt.format  = ['Absolute Loss Linear Regression|' ...
                                                'Squared Loss Linear Regression|' ...
                                                'Least Squares Regression|', ...
                                                'Student''s Loss Regression|', ...
                                                'K-Nearest Neighbour Regression|', ...
                                                'Naive Bayes Squared Loss Linear Regression|', ...
                                                'Multi-layer Perceptron with tanh Activation Function|', ...
                                                'Local Regressiony|', ...
                                                'Huberized regression|', ...
                                                'Mean regression|', ...
                                                'Stump regression|', ...
                                                'GAM model providing either smooth cubic splines, regression of degree, or linear regressions|', ...
                                                'Type-2 MLE for linear regression using Mackay''s implementation'];

                            opt.sel     = {'L1','L2','leastSquares','student','NW','NB','MLP','local','Huber','mean','stump','GAM','ARD'};

                    case 'binaryclass'
                            opt.format  = [ 'Boosting Algorithm (AdaBoost, LogitBoost)|' ...
                                            'Decision Stump Classifier|', ...
                                            'Broken Stump Classifier|', ...
                                            'Huberized hinge-loss SVM|', ...
                                            'L2-SVM|', ...
                                            'Smooth SVM with squared hinge loss|', ...
                                            'Logistic regression with L2-regularization|', ...
                                            'Elastic-net logistic regression|', ...
                                            'Random Forest|', ...
                                            'Multi-layer Perceptron|', ...
                                            'GLM with extreme link|', ...
                                            'GLM with Cauchit link|', ...
                                            'GLM optimizing the exponential loss' ];
                            opt.sel     = {'boosting','stump','brokenStump','HSVM', 'SVM', 'SSVM', 'logistic', 'elastic_net', 'randomForest', 'MLP', 'extreme', 'Cauchit','exponential'};
                end
                
                if ~isempty(param) && isfield(param,'algo')
                    f = find(strcmp(opt.sel,char(param.algo)));
                    if ~isempty(f), opt.def = f; end
                end
                
            case 'get_sublearners'
                 switch framework
                    case 'regression'
                            opt.format  = ['Naive Bayes Squared Loss Linear Regression|', ...
                                            'Mean regression|', ...
                                            'Stmup regression'];
                            opt.sel     = {'NB','mean','stump'};

                    case 'binaryclass'
                            opt.format  = [ 'Decision Stump Classifier|', ...
                                            'Broken Stump Classifier|', ...
                                            'SVM|', ...
                                            'Logistic regression with L2-regularization|', ...
                                            'GLM with extreme link|', ...
                                            'GLM with Cauchit link|', ...
                                            'GLM optimizing the exponential loss'];
                            opt.sel     = {'stump','brokenStump','SVM', 'logistic', 'extreme', 'Cauchit', 'exponential'};
                 end
                
                 if ~isempty(param) && isfield(param,'algo')
                    f = find(strcmp(opt.sel,char(param.algo)));
                    if ~isempty(f), opt.def = f; end
                 end
        end
        
    case 'get_kernel_func'
        
        switch algo
            case 'ROBSVM'
                 opt.desc    = sprintf('Select kernel function of %s', algo);
                opt.format  = 'Linear|Polynomial|RBF';
                opt.sel     = {'-t 0','-t 1','-t 2'};
                opt.def     =  1;
                opt.name    = 'kernelFunc';
            otherwise
                opt.desc    = sprintf('Select kernel function for %s', algo);
                opt.format  = 'Linear|Polynomial|RBF';
                opt.sel     = {'@ml_kernel_gram','@ml_kernel_poly','@ml_kernel_rbf'};
                opt.def     =  1;
                opt.name    = 'kernelFunc';
        end
        
    case 'get_kernel_params'

        switch algo
            case '@ml_kernel_rbf'
                opt.desc    = 'Define sigma of RBF kernel';
                opt.format  = 'e';
                opt.sel     = [];
                opt.def     = 1;
                opt.name    = 'sigma';
            case '@ml_kernel_poly'
                opt.desc    = {'Add bias to the polynomial basis', 'Define degree of polynomial basis'};
                opt.format  = {'yes|no','i'};
                opt.sel     = {[1 0],[]};
                opt.def     = {1,2};
                opt.name    = {'bias','order'};
            case '-t 1'
                opt.desc    = {'Define coefficient of polynomial basis', 'Define degree of polynomial basis'};
                opt.format  = {'i','i'};
                opt.sel     = {[],[]};
                opt.def     = {0,2};
                opt.name    = {'coef0','degree'};
            case '-t 2'
                opt.desc    = 'Define sigma of RBF kernel';
                opt.format  = 'e';
                opt.sel     = [];
                opt.def     = 2.^(-10:2:0);
                opt.name    = 'gamma';
            otherwise
                opt = [];
        end
        opt.algo = algo;
                
    case 'get_learner_params'
        
        switch algo
            case {'stump', 'brokenStump', 'mean'}
                opt = [];
            case {'leastSquares'}
                opt.desc = {'Define L2-lambda parameter range'};
                opt.format = {'e'};
                opt.sel = {[]};
                opt.def = {0};
                opt.name = {'lambdaL2'};
            case {'L1','L2','NB','SVM','exponential'}
                opt.desc = {'Add bias','Define L2-lambda parameter range'};
                opt.format = {'yes|no','e'};
                opt.sel = {[1 0],[]};
                opt.def = {1,0};
                opt.name = {'addBias','lambdaL2'};
            case {'elastic_net'}
                opt.desc = {'Add bias','Define L1-lambda parameter range', 'Define L2-lambda parameter range'};
                opt.format = {'yes|no','e','e'};
                opt.sel = {[1 0],[],[]};
                opt.def = {1,0,0};
                opt.name = {'addBias','lambdaL1','lambdaL2'};
            case 'student'
                opt.desc = {'Add bias','Define L2-lambda parameter range','Define polynomial coefficient range'};
                opt.format = {'yes|no','e','e'};
                opt.sel = {[1 0],[],[]};
                opt.def  = {1, 0, 1};
                opt.name = {'addBias','lambdaL2','poly'};
            case 'GLMNET'
                opt.desc = {'Define mixing factor (alpha) range ( 0 =ridge <-> 1 =lasso )', ...
                            'Define minimum lambda range (e.g. if N_cases>N_feats: 0.0001, otherwise 0.01)', ...
                            'Define number of lambda optimization steps', ...
                            'Define maximum number(s) of variables in the model', ...
                            'Standardize input matrix to unit variance prior to training the elastic net'};
                opt.format = {'e','e','i','i', 'yes|no'};
                opt.sel = {[],[],[],[],[1 0]};
                opt.def  = {0.5, 0.001, 100, Inf, 0};
                opt.name = {'alpha','lambda_min','nlambda','dfmax','standardize'};
            case 'GRDBST'
                opt.desc = {'Define maximum number of boosting iterations', ...
                            'Select type of loss (exponential/logarthmic for classification, squared for regression)', ...
                            'Define shrinkage factor (0<->1)', ...
                            'Define subsampling factor (0<->1)', ...
                            'Define maximum tree depth'};
                opt.format = {'i', ...
                            'exponential|logarithmic|squared', ...
                            'e', ...
                            'e', ...
                            'e'};
                opt.sel    = {[], ...
                            {'exploss','logloss','squaredloss'}, ...
                            [],...
                            [],...
                            []};
                opt.def  = {100, 1, 0.1, 0.5, 2};
                opt.name = {'maxIters','loss','shrinkageFactor','subsamplingFactor','maxTreeDepth'};
            case 'ROBSVM'
                 switch framework
                    case 'binaryclass'
                        learn_desc = 'C-SVC|nu-SVC'; learn_sel = {'-s 0', '-s 1'};

                    case 'regression'
                        learn_desc = 'epsilon-SVR|nu-SVR'; learn_sel = {'-s 3', '-s 4'}; 
                 end
                learn_def = learn_sel{1};
                opt.desc    = {'Define number of training data splits', ...
                                'Select winsorization threshold(s)', ...
                                'Select LIBSVM learner', ...
                                'Define cost parameters (only needed for C-SVC, epsilon-SVR, nu-SVR)', ...
                                'Define nu parameters (only needed for nu-SVC and nu-SVR)', ...
                                'Define epsilon parameters (only needed for epsilon-SVR)', ...
                                'Select kernel type', ...
                                'Define the kernel parameters'};
                opt.format  = {'i','i',learn_desc, 'e','e','e','kernel_func_selector','kernel_params_selector'};
                opt.sel     = {[],[],learn_sel, [], [], [], [], []};
                opt.def     = {3, 4, learn_def, 1, 0.5, 0.1, [], []};
                opt.name    = {'nsplit','wins','learner','cost','nu','epsilon','kernelFunc','kernelOptions'};
            case 'DECTRE'
                opt.desc    = {'Define parameter opimization method', ...
                                'Select type of optimization', ...
                                'Select LIBSVM learner', ...
                                'Define cost parameters (only needed for C-SVC, epsilon-SVR, nu-SVR)', ...
                                'Define nu parameters (only needed for nu-SVC and nu-SVR)', ...
                                'Define epsilon parameters (only needed for epsilon-SVR)', ...
                                'Select kernel type', ...
                                'Define the kernel parameters'};
                opt.format  = {'Cross-validation|Hyperparameter optimization','i','e','e','e'};
                opt.sel     = {{'manual','auto'},[],[], [], [] };
                opt.def     = {1, 4, 1, 0.5, 0.1 };
                opt.name    = {'nsplit','wins','learner','cost','nu','epsilon'}; 
            case 'NW'
                %uses only kernel params
                opt.def  = [];
                opt.name = [];
            case 'MLP'
                switch framework
                    case 'binaryclass'
                        opt.desc = {'Define vector of the size of each hidden layer in the multi-layer perceptron'};
                        opt.format = {'e'};
                        opt.sel = {[]};
                        opt.def  = {[3,3,3]};
                        opt.name = {'nHidden'};
                    case 'regression'
                        opt.desc = {'Define sizes of hidden layers in the multi-layer perceptron','Define active function'};
                        opt.format = {'e','s'};
                        opt.sel = {[],'tanh'};
                        opt.def  = {[3 6 9],'tanh'};
                        opt.name = {'nHidden','activFunc'};
                end    
            case 'local'
                opt.desc = {'Select submodel learner','Define submodel options',};
                opt.format = {'model_selector','options_selector','i','s'};
                opt.sel = {[], [], [], []};
                opt.def = {[], [], 1, 1};
                opt.name = {'subModel','subOptions','k','weightingFunc'};
            case 'Huber'
                opt.desc = {'Add bias','Define L2-lambda parameter range','Define epsilon parameter range'};
                opt.format = {'yes|no','e','e'};
                opt.sel = {[1 0],[],[]};
                opt.def  = {1, 0, 1};
                opt.name = {'addBias','lambdaL2','epsilon'};
            case 'ARD'
                opt.desc = {'Add bias','Define sigma2 parameter'};
                opt.format = {'yes|no','e'};
                opt.sel = {[1 0],[]};
                opt.def = {1,0};
                opt.name = {'addBias','sigma2'};
            case 'GAM'
                opt.desc = {'Define number of iterations', ...
                            'Define function used to determine f''s', ...
                            'Define highest exponent degree for polynomial subFunc', ...
                            'Define smoothness of smooth cubic splines subFunc (0 = linear regression ... 1 = unsmoothed cublic splines)', ...
                            'Define max number of iteratons for backfitting'};
                opt.format = {'e', ...
                            'smooth cubic splines|polynomial regression|linear regresion', ...
                            'e',...
                            'e',...
                            'e'};
                opt.sel = {[],{'sp1','rg','lin'},[],[],[]};
                opt.def = {500,1,2,0.8,500};
                opt.name = {'numIter','subFunc','deg','smoothing','maxIter'};
            case 'boosting'
                opt.desc = {'Select the boosting algorithm', ...
                            'Define the number of weak learners to train', ...
                            'Select the base classifier to use', ...
                            'Define the parameters of the base classifier'};
                opt.format = {'ada|logit', 'e', 'model_selector','options_selector'};
                opt.sel    = {{'ada','logit'},[],[],[]};
                opt.def    = {1,50,[],[]};
                opt.name = {'booster','nBoosts','subModel','subOptions'};
            case 'HSVM'
                opt.desc = {'Add bias', ...
                            'Define L1-lambda parameter range', ...
                            'Define L2-lambda parameter range', ...
                            'Define epsilon parameter range', ...
                            'Enable kernelization', ...
                            'Select the kernel function', ...
                            'Define the kernel parameters'};
                opt.format = {'yes|no','e','e','e','yes|no','kernel_func_selector','kernel_params_selector'};
                opt.sel = {[1 0],[],[],[], [1 0],[],[]};
                opt.def  = {1, 0, 0, 0.5, 0,[],[]};
                opt.name = {'addBias','lambdaL1','lambdaL2','epsilon','kernel','kernelFunc','kernelOptions'};
            case 'SSVM'
                opt.desc = {'Add bias', ...
                            'Define L2-lambda parameter range', ...
                            'Enable kernelization', ...
                            'Select the kernel function', ...
                            'Define the kernel parameters'};
                opt.format = {'yes|no','e','yes|no','kernel_func_selector','kernel_params_selector'};
                opt.sel = {[1 0],[],[1 0],[],[]};
                opt.def  = {1, 1, 0,[],[]};
                opt.name = {'addBias','lambdaL2','kernel','kernelFunc','kernelOptions'};
            case {'extreme','Cauchit'}
                opt.desc = {'Choose distribution (binomial for classification)', ...
                            'Define threshold'};
                opt.format = {'binomial|poisson|gamma|inverse gaussion|normal','e'};
                opt.sel = {{'binomial','poisson','gamma','inverse gaussian', 'normal'},[]};   
                opt.def  = {1, 0.5};
                opt.name = {'dist','thresh'};             
            case 'randomForest'
                opt.desc = {'Choose number of models', ...
                            'Define size of bootstrap sample (max = 1: full training sample)', ...
                            'Define minimum leaf size (the minimum number of training data in leaf nodes)', ...
                            'Choose maximum decision tree depth', ...
                            'Define number of features the each decision tree will consider (max = 1: full feature space)'};
                opt.format = {'e','e','e','e','e'};
                opt.sel = {[],[],[],[],[]};
                opt.def = {10, 0.7, 3, 8, 0.75};
                opt.name = {'nModels','nSample','minLeafSize','maxDepth','maxFeatures'};
        end
        % Overwrite defaults if pre-existing param is available
        if ~isempty(param) && isfield(param,'Params') && ~isempty(param.Params) && ~isempty(opt)
            for i=1:numel(opt.name)
                for j=1:numel(param.Params)
                    if strcmp(param.Params(j).name, opt.name{i})
                        opt.def{i} = nk_matLearn_FindDefInParam_config(opt.format{i}, opt.name{i}, opt.def{i}, opt.sel{i}, param.Params(j));
                    end
                end
            end
        end
        
end