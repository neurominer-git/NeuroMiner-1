function [X,Y,T,auc,optrocpt, indexopt, subY,subYnames] = ...
    perfcurve2(labels,scores,posClass,varargin)
%PERFCURVE Compute Receiver Operating Characteristic (ROC) curve or other
%   performance curve for classifier output.
%
%   [X,Y] = PERFCURVE(LABELS,SCORES,POSCLASS) computes a ROC curve for a
%   vector of classifier predictions SCORES given true class labels,
%   LABELS. The labels can be a numeric vector, logical vector, character
%   matrix, cell array of strings or categorical vector (see help for
%   groupingvariable). SCORES is a numeric vector of scores returned by a
%   classifier for some data. This vector must have as many elements as
%   LABELS does. POSCLASS is the positive class label (scalar), either
%   numeric (for numeric LABELS), logical (for logical LABELS) or char. The
%   specified positive class must be in the array of input labels. The
%   returned values X and Y are coordinates for the performance curve and
%   can be visualized  with PLOT(X,Y). By default, X is false positive
%   rate, FPR, (equivalently, fallout, or 1-specificity) and Y is true
%   positive rate, TPR, (equivalently, recall, or sensitivity).
%
%   [X,Y,T] = PERFCURVE(LABELS,SCORES,POSCLASS) returns an array T of
%   thresholds on classifier scores for the computed values of X and Y. It
%   has the same number of rows as X and Y. For each threshold, TP is the
%   count of true positive observations with scores greater or equal to
%   this threshold, and FP is the count of false positive observations with
%   scores greater or equal to this threshold. PERFCURVE defines negative
%   counts, TN and FN, in a similar way and sorts the thresholds in the
%   descending order which corresponds to the ascending order of positive
%   counts. For the M distinct thresholds found in the array of scores,
%   PERFCURVE returns the X, Y and T arrays with M+1 rows. PERFCURVE
%   sets elements T(2:M+1) to the distinct thresholds, and T(1) replicates
%   T(2). By convention, T(1) represents the highest 'reject all' threshold
%   and PERFCURVE computes the corresponding values of X and Y for TP=0 and
%   FP=0. T(end) is the lowest 'accept all' threshold for which TN=0 and
%   FN=0.
%
%   [X,Y,T,AUC] = PERFCURVE(LABELS,SCORES,POSCLASS) returns the area under
%   curve (AUC) for the computed values of X and Y. Unless you specify
%   XVALS or TVALS, PERFCURVE computes AUC using the returned X and Y
%   values. If XVALS or TVALS is a numeric array, PERFCURVE computes AUC
%   using X and Y values found from all distinct scores in the interval
%   specified by the smallest and largest elements of XVALS or TVALS. For
%   example, if XVALS is a numeric array, PERFCURVE finds X values for all
%   distinct thresholds and uses a subset of these (with corresponding Y
%   values) between MIN(XVALS) and MAX(XVALS) to compute AUC. The function
%   uses trapezoidal approximation to estimate the area.
%
%   If the first or last value of X or Y are NaN's, PERFCURVE removes them
%   to allow calculation of AUC. This takes care of criteria that produce
%   NaN's for the special 'reject all' or 'accept all' thresholds, for
%   example, positive predictive value (PPV) or negative predictive value
%   (NPV).
%
%   [X,Y,T,AUC] = PERFCURVE(LABELS,SCORES,POSCLASS) also returns pointwise
%   confidence bounds for the computed values X, Y, T, and AUC if you
%   supply cell arrays for LABELS and SCORES or set NBOOT to a positive
%   integer. To compute the confidence bounds, PERFCURVE uses either
%   vertical averaging (VA) or threshold averaging (TA). The returned
%   values Y are an M-by-3 array in which the 1st element in every row
%   gives the mean value, the 2nd element gives the lower bound and the 3rd
%   element gives the upper bound. The returned AUC is a row-vector with 3
%   elements following the same convention. For VA, the returned values T
%   are an M-by-3 array and X is a column-vector. For TA, the returned
%   values X are an M-by-3 matrix and T is a column-vector.
%
%   PERFCURVE computes confidence bounds using either cross-validation or
%   bootstrap. If you supply cell arrays for LABELS and SCORES, PERFCURVE
%   uses cross-validation and treats elements in the cell arrays as
%   cross-validation folds. LABELS can be a cell array of numeric vectors,
%   logical vectors, character matrices, cell arrays of strings or
%   categorical vectors. All elements in LABELS must have the same type.
%   SCORES is a cell array of numeric vectors. The cell arrays for LABELS
%   and SCORES must have the same number of elements, and the number of
%   labels in cell K must be equal to the number of scores in cell K for
%   any K in the range from 1 to the number of elements in SCORES.
%
%   If you set NBOOT to a positive integer, PERFCURVE generates NBOOT
%   bootstrap replicas to compute pointwise confidence bounds. You cannot
%   supply cell arrays for LABELS and SCORES and set NBOOT to a positive
%   integer at the same time.
%
%   PERFCURVE returns pointwise confidence bounds. It does not return
%   a simultaneous confidence band for the entire curve.
%
%   If you use 'XCrit' or 'YCrit' options described below to set the
%   criterion for X or Y to an anonymous function, PERFCURVE can only
%   compute confidence bounds by bootstrap.
%
%   [X,Y,T,AUC,OPTROCPT] = PERFCURVE(LABELS,SCORES,POSCLASS) returns the
%   optimal operating point of the ROC curve as an array of size 1-by-2
%   with FPR and TPR values for the optimal ROC operating point. OPTROCPT
%   is computed only for the standard ROC curve and set to NaN's otherwise.
%   To obtain the optimal operating point for the ROC curve, PERFCURVE
%   first finds the slope, S, using
%          S = (cost(P|N)-cost(N|N))/(cost(N|P)-cost(P|P)) * N/P
%   where cost(I|J) is the cost of assigning an observation of class J to
%   class I, and P=TP+FN and N=TN+FP are the total observation counts in
%   the positive and negative class, respectively. PERFCURVE then finds the
%   optimal operating point by moving the straight line with slope S from
%   the upper left corner of the ROC plot (FPR=0,TPR=1) down and to the
%   right until it intersects the ROC curve.
%
%   [X,Y,T,AUC,OPTROCPT,SUBY] = PERFCURVE(LABELS,SCORES,POSCLASS) returns
%   an array of Y values for negative subclasses. If you only specify one
%   negative class, SUBY is identical to Y. Otherwise SUBY is a matrix of
%   size M-by-K where M is the number of returned values for X and Y, and K
%   is the number of negative classes. PERFCURVE computes Y values by
%   summing counts over all negative classes. SUBY gives values of the Y
%   criterion for each negative class separately. For each negative class,
%   PERFCURVE places a new column in SUBY and fills it with Y values for TN
%   and FP counted just for this class.
%
%   [X,Y,T,AUC,OPTROCPT,SUBY,SUBYNAMES] = PERFCURVE(LABELS,SCORES,POSCLASS,'NegClass',NEGCLASS)
%   returns a cell array of negative class names. If you provide an input
%   array, NEGCLASS, of negative class names, PERFCURVE copies it into
%   SUBYNAMES. If you do not provide NEGCLASS, PERFCURVE extracts SUBYNAMES
%   from input labels. The order of SUBYNAMES is the same as the order of
%   columns in SUBY, that is, SUBY(:,1) is for negative class SUBYNAMES{1}
%   etc.
%
%   [X,Y] = PERFCURVE(LABELS,SCORES,POSCLASS,'PARAM1',val1,'PARAM2',val2,...)
%   specifies optional parameter name/value pairs:
%
%     'NegClass' - List of negative classes. Can be either a numeric array
%                  or an array of chars or a cell array of strings. By
%                  default, NegClass is set to 'all' and all classes found
%                  in the input array of labels that are not the positive
%                  class are considered negative. If NegClass is a subset
%                  of the classes found in the input array of labels,
%                  observations with labels that do not belong to either
%                  positive or negative classes are discarded.
%
%     'XCrit' - Criterion to compute for X. This criterion must be a
%               monotone function of the positive class score. The
%               following criteria are supported:
%         TP    - number of true positives
%         FN    - number of false negatives
%         FP    - number of false positives
%         TN    - number of true negatives
%         TP+FP - sum of TP and FP
%         RPP   = (TP+FP)/(TP+FN+FP+TN) rate of positive predictions
%         RNP   = (TN+FN)/(TP+FN+FP+TN) rate of negative predictions
%         accu  = (TP+TN)/(TP+FN+FP+TN) accuracy
%         TPR, sens, reca = TP/(TP+FN) true positive rate, sensitivity, recall
%         FNR, miss       = FN/(TP+FN) false negative rate, miss
%         FPR, fall       = FP/(TN+FP) false positive rate, fallout
%         TNR, spec       = TN/(TN+FP) true negative rate, specificity
%         PPV, prec = TP/(TP+FP) positive predictive value, precision
%         NPV       = TN/(TN+FN) negative predictive value
%         ecost=(TP*COST(P|P)+FN*COST(N|P)+FP*COST(P|N)+TN*COST(N|N))/(TP+FN+FP+TN)
%              expected cost
%         In addition, you can define an arbitrary criterion by supplying
%         an anonymous function of 3 arguments, (C,scale,cost), where C is
%         a 2-by-2 confusion matrix, scale is a 2-by-1 array of class
%         scales, and cost is a 2-by-2 misclassification cost matrix. See
%         doc for Performance Curves for more info. Warning: some of these
%         criteria return NaN values at one of the two special thresholds,
%         'reject all' and 'accept all'.
%
%     'YCrit' - Criterion to compute for Y. The same criteria as for X
%               are supported. This criterion does not have to be a
%               monotone function of the positive class score.
%
%     'XVals' - Values for the X criterion. By default, XVals is unset
%               and PERFCURVE computes X, Y and T values for all scores.
%               You can set XVals to either 'all' or a numeric array. If
%               XVals is set to 'all' and TVals is unset, PERFCURVE returns
%               X, Y and T values for all scores and computes pointwise
%               confidence bounds for Y and T using vertical averaging. If
%               XVals is set to a numeric array, PERFCURVE returns X, Y and
%               T values for the specified XVals and computes pointwise
%               confidence bounds for Y and T at these XVals using vertical
%               averaging. You cannot set XVals and TVals at the same time.
%
%     'TVals' - Thresholds for the positive class score. By default, TVals
%               is unset and PERFCURVE computes X, Y and T values for all
%               scores. You can set TVals to either 'all' or a numeric
%               array. If TVals is set to 'all' or unset and XVals is
%               unset, PERFCURVE returns X, Y and T values for all scores
%               and computes pointwise confidence bounds for Y and X using
%               threshold averaging. If TVals is set to a numeric array,
%               PERFCURVE returns X, Y and T values for the specified
%               thresholds and computes pointwise confidence bounds for Y
%               and X at these thresholds using threshold averaging. You
%               cannot set XVals and TVals at the same time.
%
%     'UseNearest' - 'on' to use nearest values found in the data instead
%                    of the specified numeric XVals or TVals and 'off'
%                    otherwise. If you specify numeric XVals and set
%                    UseNearest to 'on', PERFCURVE returns nearest unique
%                    values X found in the data, as well as corresponding
%                    values of Y and T. If you specify numeric XVals and
%                    set UseNearest to 'off', PERFCURVE returns these XVals
%                    sorted. By default this parameter is set to 'on'. If
%                    you compute confidence bounds by cross-validation or
%                    bootstrap, this parameter is always 'off'.
%
%     'ProcessNaN' - This argument specifies how PERFCURVE processes NaN
%                    scores. By default, it is set to 'ignore' and
%                    observations with NaN scores are removed from the
%                    data. If the parameter is set to 'addtofalse',
%                    PERFCURVE adds observations with NaN scores to false
%                    classification counts in the respective class. That
%                    is, observations from the positive class are always
%                    counted as false negative (FN), and observations from
%                    the negative class are always counted as false
%                    positive (FP).
%
%     'Prior' - Either string or array with 2 elements. It represents prior
%               probabilities for the positive and negative class,
%               respectively. Default is 'empirical', that is, prior
%               probabilities are derived from class frequencies. If set to
%               'uniform', all prior probabilities are set equal.
%
%     'Cost'  - A 2-by-2 matrix of misclassification costs
%                   [C(P|P) C(N|P); C(P|N) C(N|N)]
%               where C(I|J) is the cost of misclassifying
%               class J as class I. By default set to [0 0.5; 0.5 0].
%
%     'Alpha' - A numeric value between 0 and 1. PERFCURVE returns
%               100*(1-Alpha) percent pointwise confidence bounds for X, Y,
%               T and AUC. By default set to 0.05 for 95% confidence
%               interval.
%
%     'Weights' - A numeric vector of non-negative observation weights.
%                 This vector must have as many elements as SCORES or
%                 LABELS do. If you supply cell arrays for SCORES and
%                 LABELS and you need to supply WEIGHTS, you must supply
%                 them as a cell array too. In this case, every element in
%                 WEIGHTS must be a numeric vector with as many elements as
%                 the corresponding element in SCORES:
%                 NUMEL(WEIGHTS{1})==NUMEL(SCORES{1}) etc. To compute X, Y
%                 and T or to compute confidence bounds by
%                 cross-validation, PERFCURVE uses these observation
%                 weights instead of observation counts. To compute
%                 confidence bounds by bootstrap, PERFCURVE samples N out
%                 of N with replacement using these weights as multinomial
%                 sampling probabilities.
%
%     'NBoot' - Number of bootstrap replicas for computation of confidence
%               bounds. Must be a positive integer. By default this
%               parameter is set to zero, and bootstrap confidence bounds
%               are not computed. If you supply cell arrays for LABELS and
%               SCORES, this parameter must be set to zero because
%               PERFCURVE cannot use both cross-validation and bootstrap to
%               compute confidence bounds.
%
%     'BootType' - Confidence interval type used by BOOTCI to compute
%                  confidence  bounds. You can specify any type supported
%                  by BOOTCI. 'doc bootci' for more info. By default set to
%                  'bca'.
%
%     'BootArg' - Optional input arguments for BOOTCI used to compute
%                 confidence bounds. You can specify all arguments
%                 supported by BOOTCI. 'doc bootci' for more info. Empty by
%                 default.
%
%   Example: Plot ROC curve for classification by logistic regression
%      load fisheriris
%      x = meas(51:end,1:2);        % iris data, 2 classes and 2 features
%      y = (1:100)'>50;             % versicolor=0, virginica=1
%      b = glmfit(x,y,'binomial');  % logistic regression
%      p = glmval(b,x,'logit');     % get fitted probabilities for scores
%
%      [X,Y] = perfcurve(species(51:end,:),p,'virginica');
%      plot(X,Y)
%      xlabel('False positive rate'); ylabel('True positive rate')
%      title('ROC for classification by logistic regression')
%
%      % Obtain errors on TPR by vertical averaging
%      [X,Y] = perfcurve(species(51:end,:),p,'virginica','nboot',1000,'xvals','all');
%      errorbar(X,Y(:,1),Y(:,2)-Y(:,1),Y(:,3)-Y(:,1)); % plot errors
%
%   See also GLMFIT, CLASSIFY, NAIVEBAYES, CLASSREGTREE, GROUPINGVARIABLE,
%   BOOTCI.

%   Copyright 2008-2009 The MathWorks, Inc.
%   $Revision: 1.1.8.3 $

% Prepare input parser
p = inputParser;
p.addRequired('labels', ...
    @(x) ~isempty(x) && (ischar(x) || isnumeric(x) || islogical(x) ...
    || iscell(x) || strcmp(class(x),'nominal') || strcmp(class(x),'ordinal')));
p.addRequired('scores',@(x) ~isempty(x) && (isnumeric(x) || iscell(x)));
p.addRequired('PosClass', ...
    @(x) ~isempty(x) && (ischar(x) || isnumeric(x) || islogical(x)));
p.addParamValue('NegClass','all', ...
    @(x) ~isempty(x) && (ischar(x) || isnumeric(x) || iscell(x)));
p.addParamValue('XCrit','FPR',@(x) ischar(x) || isa(x,'function_handle'));
p.addParamValue('YCrit','TPR',@(x) ischar(x) || isa(x,'function_handle'));
p.addParamValue('XVals','', ...
    @(x) (ischar(x) && (strcmpi(x,'all') || isempty(x))) || ...
    (~isempty(x) && isnumeric(x)) );
p.addParamValue('TVals','', ...
    @(x) (ischar(x) && (strcmpi(x,'all') || isempty(x))) || ...
    (~isempty(x) && isnumeric(x)) );
p.addParamValue('UseNearest','on',@(x) strcmpi(x,'on') || strcmpi(x,'off'));
p.addParamValue('ProcessNaN','ignore', ...
    @(x) strcmpi(x,'ignore') || strcmpi(x,'addtofalse'));
p.addParamValue('Prior','empirical', ...
    @(x) (ischar(x) && (strcmpi(x,'empirical') || strcmpi(x,'uniform'))) ...
    || (isnumeric(x) && numel(x)==2) );
p.addParamValue('Cost',[0 0.5; 0.5 0], ...
    @(x) isnumeric(x) && size(x,1)==2 && size(x,2)==2);
p.addParamValue('Weights',[],@(x) isempty(x) || isfloat(x) || iscell(x));
p.addParamValue('NBoot',0,@(x) isnumeric(x) && isscalar(x) && x>=0);
p.addParamValue('BootType','bca',@(x) ischar(x));
p.addParamValue('BootArg',{},@(x) iscell(x));
p.addParamValue('Alpha',0.05,@(x) isnumeric(x) && isscalar(x) && x>0 && x<1);
p.FunctionName = 'perfcurve';

% Parse inputs
p.parse(labels,scores,posClass,varargin{:});
negClass       = p.Results.NegClass;
xCrit          = p.Results.XCrit;
yCrit          = p.Results.YCrit;
xVals          = p.Results.XVals;
tVals          = p.Results.TVals;
toNearest      = strcmpi(p.Results.UseNearest,'on');
processNaN     = p.Results.ProcessNaN;
prior          = p.Results.Prior;
cost           = p.Results.Cost;
weights        = p.Results.Weights;
nboot          = ceil(p.Results.NBoot);
boottype       = p.Results.BootType;
bootarg        = p.Results.BootArg;
alpha          = p.Results.Alpha;

% Prepare data
[scores,labels,weights,ncv] = preparedata(scores,labels,weights);

% Check compatibility of thresholds and xvals arguments
% By default use supplied thresholds for computing the curve
useTVals = true;
if ~isempty(tVals) && ~isempty(xVals)
    error(message('stats:perfcurve:BothTandXsupplied'));
end
% If X values are supplied, use them to compute the curve
if ~isempty(xVals)
    useTVals = false;
end

% Check if both cross-validation and bootstrap are requested
docv = false;
doboot = false;
if ncv>0
    docv = true;
end
if nboot>0
    doboot = true;
end
if docv && doboot
    error(message('stats:perfcurve:BothCVandBootstrapRequested'));
end
nsub = max(ncv,nboot);

% Stack all CV values into one long vector
% ncvel is the number of elements per CV fold
if docv
    [ncvel,scores,labels,weights] = stackcv(scores,labels,weights);
    ncvcum = [0; cumsum(ncvel)];
end

% Convert class labels to a cat array
labels = nominal(labels);
trueNames = getlabels(labels);
if length(trueNames) < 2
    error(message('stats:perfcurve:NotEnoughClasses'));
end

% Check costs
if (cost(2,1)-cost(2,2))<=0 || (cost(1,2)-cost(1,1))<=0
    error(message('stats:perfcurve:InvalidCost'));
end

% Sort scores in the descending order
[sScores,sorted] = sort(scores,1,'descend');

% Get class membership for instances:
% W(i,j) is weight of observation i if observation i is from class j and 0
% otherwise.
% Also, get negative class names.
% W has the size of NxK,
%   where N is the number of instances and K is the number of classes.
% subYnames is a cell array of length K-1 with names of negative classes.
% Column W(:,j) is for class subYnames{j-1}   (j>1)
[W,subYnames] = membership(labels(sorted),weights(sorted),...
    posClass,negClass,trueNames);

% Make Wcum, a matrix of cumulative weights in each class.
% Adjust Wcum and scores using the specified behavior for NaN scores.
% Output sorted distinct scores into sScores and corresponding rows into Wcum.
% Wcum and output sScores do not have the same size as W and input sScores.
% To access the full vector of scores, use scores(sorted).
[Wcum,sScores] = makeccum(W,sScores,processNaN);

% Get cumulative counts and sorted scores for each subset
if docv
    subscores = cell(nsub,1);
    Wcumsub = cell(nsub,1);
    for isub=1:nsub
        tf = sorted>ncvcum(isub) & sorted<=ncvcum(isub+1);
        [Wcumsub{isub},subscores{isub}] = ...
            makeccum(W(tf,:),scores(sorted(tf)),processNaN);
    end
end

% Check that confidence bound computation by CV is not requested for
% user-defined criteria
if (isa(xCrit,'function_handle') || isa(yCrit,'function_handle')) && docv
    error(message('stats:perfcurve:UserCritConfBounds'));
end

% Determine criteria to compute
% The 1st output is the function for computing criterion itself.
% The 2nd output is the function for computing weight for this criterion.
[fx,fwx] = makeCrit(xCrit);
[fy,fwy] = makeCrit(yCrit);

% Convert class probabilities into class scales.
scale = classscale(Wcum,prior);

% Define arrayfuns
afx = @(tp,fn,fp,tn) fx([tp fn; fp tn],scale,cost);
afy = @(tp,fn,fp,tn) fy([tp fn; fp tn],scale,cost);
afwx = @(tp,fn,fp,tn) fwx([tp fn; fp tn],scale);
afwy = @(tp,fn,fp,tn) fwy([tp fn; fp tn],scale);

% fix for faster performance:
if strcmp(xCrit,'FPR')
    afx = @(tp,fn,fp,tn)(fp)/fp(end);
    afwx = @(tp,fn,fp,tn)fp(end)*ones(size(fp));
end
if strcmp(yCrit,'TPR')
    afy = @(tp,fn,fp,tn)(tp)/tp(end);
    afwy = @(tp,fn,fp,tn)tp(end) * ones(size(tp));

end
% Compute threshold indices
uniqdiv = ~docv && ~doboot && toNearest;
Ndiv = size(Wcum,1) - 1;
if ischar(xVals) && ischar(tVals)
    div = (1:Ndiv)'; % Use all thresholds
else
    if useTVals % Use user-defined thresholds
        [div,tVals] = tdiv(tVals,sScores,uniqdiv);
    else % Use user-defined X values
        [div,xVals] = xdiv(xVals,afx,Wcum,uniqdiv);
    end
end

% Compute the actual values for the specified criterion,
%   (to be plotted on X axis),
%   and associated TP and FP counts.
[X,~,tpX,fpX] = Xvalues(div,Wcum,afx,[]);
% If xVals are supplied by the user and if they do not need to be set to
% nearest found values, use these xVals
if ~uniqdiv && isnumeric(xVals)
    X = xVals(:);
end

% Compute criterion values associated with these thresholds
%   (to be plotted on Y axis)
Y = Yvalues(tpX,fpX,Wcum,afy,[]);

% Check if any criteria are NaN's besides those at special 'reject all' and
% 'accept all' thresholds
special = (div==1 | div==Ndiv);
if any(isnan(X(~special)))
    error(message('stats:perfcurve:BadXCritValue'));
end
if any(isnan(Y(~special)))
    error(message('stats:perfcurve:BadYCritValue'));
end

% Get thresholds from indices
if nargout>2
    % If tVals are supplied by the user and if they do not need to be set
    % to nearest found values, use these tVals
    if ~uniqdiv && isnumeric(tVals)
        T = tVals(:);
    else
        T = thresholds(div,sScores,false);
    end
end

% Find numeric divisions for VA or TA
tValsFixed = [];
xValsFixed = [];
if docv || doboot
    if useTVals
        if ischar(tVals)
            tValsFixed = thresholds(div,sScores,true);
        else
            tValsFixed = tVals;
        end
    else
        if ischar(xVals)
            xValsFixed = X;
        else
            xValsFixed = xVals;
        end
    end
end

% Compute confidence intervals
ciX = [];
ciY = [];
ciT = [];
if docv || doboot
    notifyClassAbsent(true); % Reset notification flag.
    if     docv % cross-validation
        % For CV, compute X, Y and T for supplied folds with weights
        [Xsub,wXsub,Ysub,wYsub,Tsub,wTsub] = ...
            xyt(xValsFixed,tValsFixed,subscores,Wcumsub,afx,afwx,afy,afwy);
        
        % X errors for TA
        if useTVals
            ciX = cvci(X,Xsub,wXsub,alpha);
        end
        
        % Score threshold errors for VA
        computeT = nargout>2 && ~useTVals;
        if computeT
            ciT = cvci(T,Tsub,wTsub,alpha);
        end
        
        % Always need Y errors
        ciY = cvci(Y,Ysub,wYsub,alpha);
    elseif doboot
        % Use bootci to compute confidence intervals. Use weights as
        % multinomial probabilities for bootstrap replica generation and
        % use observation counts (not observation weights!) to compute
        % performance curves on these replicas.
        bootfun = @(idx) ...
            oneBootXYT(idx,scores(sorted),W>0,xValsFixed,tValsFixed,processNaN,afx,afy);
        ci = bootci(nboot,{bootfun (1:size(W,1))'},...
            'weights',sum(W,2),'alpha',alpha,'type',boottype,bootarg{:});
        
        % X and Y errors for TA
        if useTVals
            if length(tValsFixed)>1
                ciX = [ci(1,:,1)' ci(2,:,1)'];
                ciY = [ci(1,:,2)' ci(2,:,2)'];
            else
                ciX = ci(:,1)';
                ciY = ci(:,2)';
            end
            % Y errors for VA
        else
            if length(xValsFixed)>1
                ciY = [ci(1,:,2)' ci(2,:,2)'];
            else
                ciY = ci(:,2)';
            end
        end
        
        % Score threshold errors for VA
        computeT = nargout>2 && ~useTVals;
        if computeT
            if length(xValsFixed)>1
                ciT = [ci(1,:,1)' ci(2,:,1)'];
            else
                ciT = ci(:,1)';
            end
        end
    end
end

% Compute area under curve
if nargout>3
    % Set range and either 'x' or 't' to specify what range
    XorTrange = [];
    mode = '';
    if     ~ischar(xVals)
        mode = 'x';
        XorTrange = sort([xVals(1) xVals(end)]);
    elseif ~ischar(tVals)
        mode = 't';
        XorTrange = sort([tVals(1) tVals(end)]);
    end
    
    % Get area
    auc = AUC(XorTrange,mode,sScores,Wcum,afx,afy);
    
    % Get error
    if docv || doboot
        if     docv
            subauc = zeros(1,nsub);
            Wtotsub = zeros(1,nsub);
            for isub=1:nsub
                [subauc(isub),Wtotsub(isub)] = ...
                    AUC(XorTrange,mode,subscores{isub},Wcumsub{isub},afx,afy);
            end
            ciauc = cvci(auc,subauc,Wtotsub,alpha);
        elseif doboot
            bootfun = @(idx) ...
                oneBootAUC(idx,scores(sorted),W>0,XorTrange,mode,processNaN,afx,afy);
            ciauc = bootci(nboot,{bootfun (1:size(W,1))'},...
                'weights',sum(W,2),'alpha',alpha,'type',boottype,bootarg{:});
            ciauc = ciauc';
        end
        auc = [auc ciauc];
    end
end

% Find optimal operation point for the standard ROC curve
if nargout>4
    isroc = (strcmpi(xCrit,'FPR') || strcmpi(xCrit,'fall')) && ...
        (strcmpi(yCrit,'TPR') || strcmpi(yCrit,'sens') || strcmpi(yCrit,'reca'));
    if isroc
        [optrocpt,indexopt] = findoptroc(X,Y,Wcum,scale,cost);
    else % Not a standard ROC curve.
        optrocpt = NaN(1,2); indexopt = NaN;
    end
end

% Compute criterion values for individual negative classes
if nargout>5
    subY = subYvalues(tpX,fpX,div,Wcum,afy);
end

% Include confidence intervals if they were computed
if ~isempty(ciX)
    X = [X ciX];
end
if ~isempty(ciY)
    Y = [Y ciY];
end
if ~isempty(ciT)
    T = [T ciT];
end
end


function [W,negClassNames] = membership(sLabels,sWeights,posClass,negClass,trueNames)
% Find the positive class. Must have exactly one.
posClass = cellstr(nominal(posClass));
if length(posClass)>1
    error(message('stats:perfcurve:TooManyPositiveClasses'));
end
if ~ismember(posClass,trueNames)
    error(message('stats:perfcurve:PositiveClassNotFound'));
end

% Check negative class labels
if strcmpi(negClass,'all')
    negClass = nominal(trueNames);
    [~,posClassLoc] = ismember(posClass,cellstr(negClass));
    negClass(posClassLoc) = [];
    negClass = nominal(negClass);
else
    negClass = nominal(negClass);
    tf = ismember(cellstr(negClass),trueNames);
    if any(~tf)
        error(message('stats:perfcurve:NegativeClassNotFound'));
    end
    tf = ismember(posClass,cellstr(negClass));
    if tf
        error(message('stats:perfcurve:PositiveAndNegativeClassesOverlap'));
    end
end
nNeg = length(negClass);

% Names of selected negative classes
negClassNames = getlabels(negClass);

% Check for duplicate entries
if nNeg~=length(negClassNames)
    error(message('stats:perfcurve:DuplicateNegativeClasses'));
end

% Fill out the membership matrix
% 1st column is for the positive class.
% Columns 2:end are for negative classes.
negClass = cellstr(negClass);
C = false(length(sLabels),1+nNeg);
C(:,1) = ismember(sLabels,posClass);
for i=1:nNeg
    C(:,i+1) = ismember(sLabels,negClass(i));
end

% Get weighted membership matrix
W = bsxfun(@times,C,sWeights);
end


function [Wcum,scores] = makeccum(W,scores,processNaN)
% Discard instances that do not belong to any class
idxNone = ~any(W,2);
W(idxNone,:) = [];
scores(idxNone) = [];

% Get rid of NaN's in scores
Wnanrow = zeros(1,size(W,2));
idxNaN = isnan(scores);
if strcmpi(processNaN,'addtofalse')
    if ~isempty(idxNaN)
        Wnanrow = sum(W(idxNaN,:),1);
    end
end
scores(idxNaN) = [];
W(idxNaN,:) = [];

% Make a matrix of counts with NaN instances included
Wnan = zeros(size(W,1)+2,size(W,2));
Wnan(1,2:end) = Wnanrow(2:end);% FP (always accepted)
Wnan(2:end-1,:) = W;
Wnan(end,1) = Wnanrow(1);% FN (always rejected)

% Compute cumulative counts in each class
Wcum = cumsum(Wnan,1);

% Compact Wcum in case of identical scores
idxEq = find( scores(1:end-1) < scores(2:end) + ...
    max([eps(scores(1:end-1)) eps(scores(2:end))],[],2) );
Wcum(idxEq+1,:) = [];
scores(idxEq) = [];
end


function scale = classscale(Wcum,prior)
scale = zeros(2,1);
Wpos = Wcum(end,1);
Wneg = sum(Wcum(end,2:end),2);
if ischar(prior) && strcmpi(prior,'empirical')
    scale = ones(2,1);
end
if ischar(prior) && strcmpi(prior,'uniform')
    prior = ones(2,1);
end
if isnumeric(prior)
    if any(prior<=0)
        error(message('stats:perfcurve:NonPositivePriors'));
    end
    scale(1) = prior(1)*Wneg;
    scale(2) = prior(2)*Wpos;
    scale = scale/sum(scale);
end
end


function [f,wf] = makeCrit(crit)
% If this is a user-supplied function, just return it
if isa(crit,'function_handle')
    f = crit;
    wf = @(C,scale) sum(scale.*sum(C,2));
    return;
end

% Make the function given criterion name
switch lower(crit)
    case 'tp'
        f = @(C,scale,cost) scale(1)*C(1,1);
        wf = @(C,scale) 1;
    case 'fn'
        f = @(C,scale,cost) scale(1)*C(1,2);
        wf = @(C,scale) 1;
    case 'fp'
        f = @(C,scale,cost) scale(2)*C(2,1);
        wf = @(C,scale) 1;
    case 'tn'
        f = @(C,scale,cost) scale(2)*C(2,2);
        wf = @(C,scale) 1;
    case 'tp+fp'
        f = @(C,scale,cost) sum(scale.*C(:,1),1);
        wf = @(C,scale) 1;
    case 'rpp'
        f = @(C,scale,cost) sum(scale.*C(:,1),1) / sum(scale.*sum(C,2));
        wf = @(C,scale) sum(scale.*sum(C,2));
    case 'rnp'
        f = @(C,scale,cost) sum(scale.*C(:,2),1) / sum(scale.*sum(C,2));
        wf = @(C,scale) sum(scale.*sum(C,2));
    case 'accu'
        f = @(C,scale,cost) ...
            (scale(1)*C(1,1)+scale(2)*C(2,2)) / sum(scale.*sum(C,2));
        wf = @(C,scale) sum(scale.*sum(C,2));
    case {'tpr','sens','reca'}
        f = @(C,scale,cost) C(1,1) / sum(C(1,:),2);
        wf = @(C,scale) sum(C(1,:),2);
    case {'fnr','miss'}
        f = @(C,scale,cost) C(1,2) / sum(C(1,:),2);
        wf = @(C,scale) sum(C(1,:),2);
    case {'fpr','fall'}
        f = @(C,scale,cost) C(2,1) / sum(C(2,:),2);
        wf = @(C,scale) sum(C(2,:),2);
    case {'tnr','spec'}
        f = @(C,scale,cost) C(2,2) / sum(C(2,:),2);
        wf = @(C,scale) sum(C(2,:),2);
    case {'ppv','prec'}
        f = @(C,scale,cost) scale(1)*C(1,1) / sum(scale.*C(:,1));
        wf = @(C,scale) sum(scale.*C(:,1));
    case 'npv'
        f = @(C,scale,cost) scale(2)*C(2,2) / sum(scale.*C(:,2));
        wf = @(C,scale) sum(scale.*C(:,2));
    case 'ecost'
        f = @(C,scale,cost) ...
            sum(scale.*sum(C.*cost,2)) / sum(scale.*sum(C,2));
        wf = @(C,scale) sum(scale.*sum(C,2));
    otherwise
        error(message('stats:perfcurve:UnknownXYCriterion'));
end
end


function inrange = applyrange(X,T,XorTrange,mode)
inrange = true(length(X),1);
if ~isempty(XorTrange)
    if numel(XorTrange)~=2
        error(message('stats:perfcurve:InvalidXorTRange'));
    end
    XorTrange = sort(XorTrange);
    if     mode=='x'
        inrange = (X>=XorTrange(1) & X<=XorTrange(2));
    elseif mode=='t'
        inrange = (T>=XorTrange(1) & T<=XorTrange(2));
    end
    if isempty(inrange)
        error(message('stats:perfcurve:XorTRangeTooRestrictive'));
    end
end
end


function T = thresholds(div,scores,shiftRejectAll)
% Init
T = zeros(length(div),1);

% First threshold
isone = div==1;
T(isone) = scores(1);
if shiftRejectAll
    T(isone) = T(isone) + eps(scores(1));
end

% Normal thresholds
T(~isone) = scores(div(~isone)-1);
end


function increasing = monotone(vals)
% Allow only monotone criteria.
% By default, assume a criterion that monotonously increases as
%   the predicted score in the positive class decreases. Otherwise, swap.
% 'increasing' is a flag that shows in what order values are sorted.
vals = vals(~isnan(vals));
increasing = 1;
if isempty(vals)
    return;
end
if any(vals(1:end-1)>vals(2:end))
    increasing = -1;
    if any(vals(1:end-1)<vals(2:end))
        error(message('stats:perfcurve:NonMonotoneXCriterion'));
    end
end
end

%{
% This routine does the same thing as the uncommented one below. It is
% simpler but could be much slower for huge allVals.
function div = finddiv(divVals,allVals,increasing)
Nall = length(allVals);
Nval = length(divVals);
div = ones(Nval,1);
for i=1:Nall
    div(increasing*divVals >= increasing*allVals(i)) = i;
end
end
%}

% Find allVals indices for values supplied in divVals.
% allVals is a sorted array of all know values, and divVals is a sorted
% array of values for which we want to find indices in allVals.
% increasing is the sort order: +1 for ascending and -1 for descending.
% Assume that divVals is not too large and allVals can be huge.
% This routine works if allVals is a subset of divVals as well (case
% relevant for cross-validation or bootstrap).
function div = finddiv(divVals,allVals,increasing)
% Init
Ndiv = length(divVals);
div = ones(Ndiv,1);

% Set current search starting index
iAll = 1;

% Main loop
for iDiv=1:Ndiv
    % Find index in the list of all known values
    thisAll = find(increasing*allVals(iAll:end) <= increasing*divVals(iDiv),1,'last');
    if isempty(thisAll)
        continue;
    end
    
    % Update starting index for the list of all known values
    iAll = iAll + thisAll - 1;
    
    % Fill found divisions
    div(iDiv:end) = iAll;
end
end


% Find division indices along X for supplied xVals
function [divX,xVals] = xdiv(xVals,afx,Wcum,uniqdiv)
% Get counts for positive and negative classes
Nrow = size(Wcum,1) - 1;
Pcum = Wcum(:,1);
Ncum = sum(Wcum(:,2:end),2);
wP = Pcum(end);
wN = Ncum(end);

% Get all possible values of the criterion
try
    allVals = afx(Pcum(1:Nrow),wP-Pcum(1:Nrow),Ncum(1:Nrow),wN-Ncum(1:Nrow));
catch
    allVals = arrayfun(afx,Pcum(1:Nrow),wP-Pcum(1:Nrow),Ncum(1:Nrow),wN-Ncum(1:Nrow));
end
% Do criterion values increase or decrease vs predicted scores?
increasing = monotone(allVals);

% Sort input values
xVals = increasing*sort(increasing*xVals);

% Find indices of thresholds.
divX = finddiv(xVals,allVals,increasing);

% Make indices unique if necessary
% and include the first threshold below 'accept all'.
if uniqdiv
    if increasing*allVals(1)>=increasing*xVals(1)
        divX = [1; divX];
    end
    divX = unique(divX);
end
end


% Find division indices along positive class score for supplied thresholds
function [divT,tVals] = tdiv(tVals,scores,uniqdiv)
% Sort thresholds
tVals = sort(tVals,'descend');

% Find indices of thresholds.
% Scores have been sorted in the descending order.
divT = finddiv(tVals,scores,-1);

% Offset indices of non-special thresholds by 1 to account for the NaN row
% on top of Ccum.
divT = divT + 1;
divT(tVals>scores(1)) = 1;

% Make indices unique if necessary
if uniqdiv
    divT = unique(divT);
end
end


% Get X values for given division indices
function [valX,wX,tpX,fpX] = Xvalues(div,Wcum,afx,afwx)
% Get counts for positive and negative classes
Pcum = Wcum(:,1);
Ncum = sum(Wcum(:,2:end),2);
wP = Pcum(end);
wN = Ncum(end);
Nrow = size(Pcum,1)-1;

% Check that division index does not exceed array size
if any(div>Nrow)
    error(message('stats:perfcurve:BadDivisionIndex'));
end

% Get TP and FP counts for chosen threshold indices
tpX = Pcum(div);
fpX = Ncum(div);

% valX is the corresponding array of criterion values
try
    valX = afx(tpX,wP-tpX,fpX,wN-fpX);
catch
    valX = arrayfun(afx,tpX,wP-tpX,fpX,wN-fpX);
end

% wX is the corresponding array of weights
wX = [];
if ~isempty(afwx)
    try 
        wX = afwx(tpX,wP-tpX,fpX,wN-fpX);
    catch
        wX = arrayfun(afwx,tpX,wP-tpX,fpX,wN-fpX);
    end
end
end


function [valY,wY] = Yvalues(tpX,fpX,Wcum,afy,afwy)
% Get number of instances in the positive and negative classes
wP = Wcum(end,1);
wN = sum(Wcum(end,2:end),2);

% Compute Y criterion for the total of all negative classes
try
    valY = afy(tpX,wP-tpX,fpX,wN-fpX);
catch
    valY = arrayfun(afy,tpX,wP-tpX,fpX,wN-fpX);
end


% Compute corresponding weights
wY = [];
if ~isempty(afwy)
    try 
        wY = afwy(tpX,wP-tpX,fpX,wN-fpX);
    catch
        wY = arrayfun(afyx,tpX,wP-tpX,fpX,wN-fpX);
    end
end
end


function subY = subYvalues(tpX,fpX,divX,Wcum,afy)
% If only one negative class, it is accounted for by Yvalues
nNegClass = size(Wcum,2)-1;
if nNegClass < 2
    subY = Yvalues(tpX,fpX,Wcum,afy,[]);
    return;
end

% Compute Y criteria for negative classes separately
wP = Wcum(end,1);
subY = zeros(length(tpX),nNegClass);
for i=1:nNegClass
    wN = Wcum(end,i+1);
    if wN==0
        error(message('stats:perfcurve:BadClassCounts', i));
    end
    fpX = Wcum(divX,i+1);
    subY(:,i) = arrayfun(afy,tpX,wP-tpX,fpX,wN-fpX);
end
end


% Compute area under curve and associated weight
function [auc,wauc] = AUC(XorTrange,mode,sScores,Wcum,afx,afy)
% If not all thresholds have been found, find all and apply range.
Ndiv = size(Wcum,1)-1;
div = (1:Ndiv)';
[X,~,tpx,fpx] = Xvalues(div,Wcum,afx,[]);
Y = Yvalues(tpx,fpx,Wcum,afy,[]);
T = thresholds(div,sScores,false);

% Apply range
inrange = applyrange(X,T,XorTrange,mode);
needTrimFirst = true;
if ~inrange(1)
    needTrimFirst = false;
end
needTrimLast = true;
if ~inrange(end)
    needTrimLast = false;
end
X = X(inrange);
Y = Y(inrange);
icum = find(inrange,1,'last');
if isempty(icum)
    wauc = 0;
else
    wauc = sum(Wcum(icum,:),2);
end

% If the 1st or last value is NaN, trim the new X and Y
if needTrimLast && (isnan(X(end)) || isnan(Y(end)))
    X(end) = [];
    Y(end) = [];
end
if needTrimFirst && (isnan(X(1)) || isnan(Y(1)))
    X(1) = [];
    Y(1) = [];
end

% Have enough data?
if length(X)<2
    auc = 0;
    return;
end

% Get area
auc = 0.5*sum( (X(2:end)-X(1:end-1)).*(Y(2:end)+Y(1:end-1)) );
auc = abs(auc);
end


function [optpt, idx] = findoptroc(X,Y,Wcum,scale,cost)
% Get positive and negative counts
wP = scale(1)*Wcum(end,1);
wN = scale(2)*sum(Wcum(end,2:end),2);

% Get the optimal slope
m = (cost(2,1)-cost(2,2))/(cost(1,2)-cost(1,1)) * wN/wP;

% Find lowest intercept for straight lines drawn through (X,Y)
%   using this slope and X axis
[~,idx] = min(X - Y/m);

% Get the optimal point
optpt = [X(idx) Y(idx)];
end


function [scores,labels,weights,ncv] = preparedata(scores,labels,weights)
ncv = 0; % no cross-validation
if iscell(scores) % Cross-validated scores and labels
    % Scores
    if ~isvector(scores)
        error(message('stats:perfcurve:CVScoresNotVector'));
    end
    scores = scores(:);
    ncv = numel(scores);
    if ncv<2
        error(message('stats:perfcurve:CVScoresTooShort'));
    end
    tfisemp = cellfun(@isempty,scores);
    tfisnum = cellfun(@isnumeric,scores);
    tfisvec = cellfun(@isvector,scores);
    if any(tfisemp) || any(~tfisnum) || any(~tfisvec)
        error(message('stats:perfcurve:CVScoresWithBadElements'));
    end
    nels = cellfun(@numel,scores);
    for icv=1:ncv
        scores{icv} = scores{icv}(:);
    end
    
    % Labels
    if ~iscell(labels)
        error(message('stats:perfcurve:CVLabelsNotCell'));
    end
    if ~isvector(labels)
        error(message('stats:perfcurve:CVLabelsNotVector'));
    end
    if numel(labels)~=ncv
        error(message('stats:perfcurve:CVLabelsWithNonmatchingLength'));
    end
    labels = labels(:);
    ltype = cellfun(@class,labels,'UniformOutput',false);
    for icv=2:ncv
        if ~strcmp(ltype{1},ltype{icv})
            error(message('stats:perfcurve:CVLabelsWithDifferentTypes'));
        end
    end
    nell = zeros(ncv,1);
    for icv=1:ncv
        if strcmp(ltype{icv},'char')
            nell(icv) = size(labels{icv},1);
        else
            labels{icv} = labels{icv}(:);
            nell(icv) = numel(labels{icv});
        end
    end
    
    % Weights
    if isempty(weights)
        weights = cell(ncv,1);
    else
        if ~iscell(weights)
            error(message('stats:perfcurve:CVWeightsNotCell'));
        end
        if ~isvector(weights)
            error(message('stats:perfcurve:CVWeightsNotVector'));
        end
        if numel(weights)~=ncv
            error(message('stats:perfcurve:CVWeightsWithNonmatchingLength'));
        end
        weights = weights(:);
    end
    for icv=1:ncv
        if isempty(weights{icv})
            weights{icv} = ones(nels(icv),1);
        else
            weights{icv} = weights{icv}(:);
        end
    end
    tfisflo = cellfun(@isfloat,weights);
    tfisvec = cellfun(@isvector,weights);
    if any(~tfisflo) || any(~tfisvec)
        error(message('stats:perfcurve:CVWeightsNotNumeric'));
    end
    fneg = @(x) any(x<0);
    fzer = @(x) all(x==0);
    fnan = @(x) all(isnan(x));
    tfisneg = cellfun(fneg,weights);
    tfiszer = cellfun(fzer,weights);
    tfisnan = cellfun(fnan,weights);
    if any(tfisneg) || any(tfiszer) || any(tfisnan)
        error(message('stats:perfcurve:CVWeightsNegativeOrNaN'));
    end
    nelw = cellfun(@numel,weights);
    
    % Check sizes
    if any(nels~=nell) || any(nelw~=nell)
        error(message('stats:perfcurve:CVWeightsNotMatchedToScores'));
    end
    
    % Get rid of observations with zero weights
    for icv=1:ncv
        iszer = weights{icv}==0;
        scores{icv}(iszer) = [];
        labels{icv}(iszer,:) = [];
        weights{icv}(iszer) = [];
    end
else % Plain scores, labels and weights
    % Scores
    if ~isvector(scores)
        error(message('stats:perfcurve:ScoresNotVector'));
    end
    scores = scores(:);
    ns = numel(scores);
    
    % Labels
    if ~ischar(labels)
        if ~isvector(labels)
            error(message('stats:perfcurve:LabelsNotVector'));
        end
        labels = labels(:);
    end
    nl = size(labels,1);
    if ns~=nl
        error(message('stats:perfcurve:ScoresAndLabelsDoNotMatch'));
    end
    
    % Weights
    if isempty(weights)
        weights = ones(ns,1);
    else
        if ~isfloat(weights)
            error(message('stats:perfcurve:WeightsNotRealNumbers'));
        end
        if ~isvector(weights)
            error(message('stats:perfcurve:WeightsNotVector'));
        end
        if any(weights<0) || all(weights==0)
            error(message('stats:perfcurve:WeightsNegative'));
        end
        if numel(weights)~=ns
            error(message('stats:perfcurve:WeightsNotMatchedToScores'));
        end
        if any(isnan(weights))
            error(message('stats:perfcurve:NaNWeights'));
        end
        weights = weights(:);
        iszer = weights==0;
        scores(iszer) = [];
        labels(iszer,:) = [];
        weights(iszer) = [];
    end
end
end


% Stack cell arrays of scores used for cross-validation, labels and weights
% into long vectors preserving type. ncvel is a vector storing the number
% of elements in every CV piece.
function [ncvel,scores,labels,weights] = stackcv(scores,labels,weights)
ncvel = cellfun(@numel,scores);
scores = vertcat(scores{:});
weights = vertcat(weights{:});
if ischar(labels{1})
    labels = char(labels{:});
else
    labels = vertcat(labels{:});
end
end


% Compute division indices, and X, Y and T values for either
% cross-validation or bootstrap samples. scores and Wcum are cell arrays
% with Nsub elements, one element per subset. wX and wY are weights
% appropriate for computation of statistics X and Y. For thresholds T,
% always use the total weight of the subset for any threshold.
function [X,wX,Y,wY,T,wT] = xyt(xVals,tVals,scores,Wcum,afx,afwx,afy,afwy)
% Init
if     ~isempty(xVals)
    M = length(xVals);
else
    %elseif ~isempty(tVals)
    M = length(tVals);
end
nsub = length(scores);
X = zeros(M,nsub);
Y = zeros(M,nsub);
wX = zeros(M,nsub);
wY = zeros(M,nsub);

% Loop through subsets
T = [];
wT = [];
computeT = nargout>4 && ~isempty(xVals);
if computeT
    T = zeros(M,nsub);
    wT = zeros(M,nsub);
    for isub=1:nsub
        [X(:,isub),wX(:,isub),Y(:,isub),wY(:,isub),T(:,isub),wT(:,isub)] = ...
            xytOneSample(xVals,tVals,scores{isub},Wcum{isub},afx,afwx,afy,afwy);
    end
else
    for isub=1:nsub
        [X(:,isub),wX(:,isub),Y(:,isub),wY(:,isub)] = ...
            xytOneSample(xVals,tVals,scores{isub},Wcum{isub},afx,afwx,afy,afwy);
    end
end
end


% Returns one sample of XYT values with corresponding weights.
% Weights for T averaging are set to the cumulative weight of the sample.
function [X,wX,Y,wY,T,wT] = xytOneSample(xVals,tVals,scores,Wcum,afx,afwx,afy,afwy)
% Notify the user once that the class is missing in one of the folds.
if ~all(sum(Wcum,1)) && notifyClassAbsent()
    warning(message('stats:perfcurve:SubSampleWithMissingClasses'));
end

% Get divisions
if     ~isempty(tVals) % perform threshold averaging
    div = tdiv(tVals,scores,false);
else
    %elseif ~isempty(xVals) % perform vertical averaging
    div = xdiv(xVals,afx,Wcum,false);
end

% Get X and Y
[X,wX,tpX,fpX] = Xvalues(div,Wcum,afx,afwx);
[Y,wY] = Yvalues(tpX,fpX,Wcum,afy,afwy);

% Get T
T = [];
wT = [];
computeT = nargout>4 && ~isempty(xVals);
if computeT
    if     ~isempty(xVals)
        M = length(xVals);
    elseif ~isempty(tVals)
        M = length(tVals);
    end
    T = thresholds(div,scores,false);
    wT = repmat(sum(Wcum(end,:),2),M,1);
end
end


% Returns CV confidence interval for every row of Xcv.
function ci = cvci(X,Xcv,wXcv,alpha)
ci = [];
ncv = size(Xcv,2);
if ncv>0
    e = wstd(X,Xcv,wXcv)/sqrt(ncv);
    ea = e*norminv(1-0.5*alpha,0,1);
    ci = [X-ea X+ea];
end
end


% Returns one bootstrap replica specified by indices idx.
% Weights are not needed here because they are used to generate replicas
% and must be fed into bootci directly.
function out = oneBootXYT(idx,scores,C,xVals,tVals,processNaN,afx,afy)
idx = sort(idx);
[Ccum,subscores] = makeccum(C(idx,:),scores(idx),processNaN);
if     isempty(xVals) % TA
    [X,~,Y] = xytOneSample(xVals,tVals,subscores,Ccum,afx,[],afy,[]);
    out = [X Y];
else
    %elseif isempty(tVals) % VA
    [~,~,Y,~,T,~] = xytOneSample(xVals,tVals,subscores,Ccum,afx,[],afy,[]);
    out = [T Y];
end
end


function out = oneBootAUC(idx,scores,C,XorTrange,mode,processNaN,afx,afy)
idx = sort(idx);
[Ccum,subscores] = makeccum(C(idx,:),scores(idx),processNaN);
out = AUC(XorTrange,mode,subscores,Ccum,afx,afy);
end


% Compute weighted standard deviation given mean M,
% measurements X and their weights W.
% X and W are NxK matrices for N measurements and K folds.
function E = wstd(M,X,W)
tfnan = isnan(X);
X(tfnan) = 0;
W(tfnan) = 0;
tfzeroW = ~any(W,2);
if any(tfzeroW)
    warning(message('stats:perfcurve:CVfoldsWithZeroWeights'));
end
E = zeros(length(M),1);
i = ~tfzeroW;
E(i) = sqrt(sum(W(i,:).*bsxfun(@minus,X(i,:),M(i)).^2,2) ./ sum(W(i,:),2));
end


% Notify the user if a class is absent?
function tf = notifyClassAbsent(reset)
persistent notified;
if nargin>0 && reset
    tf = true;
    notified = false;
else
    tf = ~notified;
    notified = true;
end
end
