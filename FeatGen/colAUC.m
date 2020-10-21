function [auc, labels] = colAUC(X, y, varargin)
%COLAUC calculates Area under ROC curve (AUC) for a vector or for each
% column of a matrix.
%
% SYNTAX:
%   auc = colAUC(X, y)      - calculates AUC for each column of X
%   auc = colAUC(X, y, alg) - calculates AUC for each column of X using
%                             chosen algorithm
%   [auc labels] = colAUC(X, y) - if "y" has more than 2 clases than AUC
%             will be calculated for each pair of classes. Each row in
%             "labels" array matches a row in "auc" array and returns
%             labels of 2 classes compared, in the same format as input "y"
%   colAUC(X, y) - if no output is requested than colAUC plots ROC curves
%   auc = colAUC(X, y, ...) - takes additional named arguments:
%   'algorithm', 'plot' and 'abs'
%
% INPUTS:
%   X - 1D vector or 2D matrix (of feature columns and samples rows))
%   y - Class labels for the X data samples in a form of a vector with one
%       label for each row/sample of X. Can be either a cell array or
%       numeric vector.
%   alg - Algorithm to use: "ROC" integrates ROC curves, while "Wilcoxon"
%       uses Wilcoxon Rank Sum Test to get the same results. Default
%       "Wilcoxon" is faster. This argument is mostly provided for verification.
%   plot - boolean - plot ROC curves? By default false if function returns
%       anything and true if it doesn't.
%   abs - boolean - By default AUC returned is always between 0.5 and 1,
%       and if calculated AUC is less than 0.5 than AUC=1-AUC is applied. 
%       behavior is useful if you are using AUC to calculate power of each
%       This feature to distinguish between 2 classes without assigning 
%       meaning to each class. If that behavior is unwanted than setting 
%       'abs' to false will disable it. Than range of possible values of
%       AUC will be between 0 and 1.
%
% DISCUSION
%  AUC is a very useful measure of similarity between two classes
%  measuring area under "Receiver Operating Characteristic" or ROC curve.
%  In case of data with no ties all sections of ROC curve are either
%  horizontal or vertical, in case of data with ties diagonal
%  sections can also occur. Area under the ROC curve is calculated using
%  trapz function. AUC is always in between 0.5 (two classes are
%  statistically identical) and 1.0 (there is a threshold value that can
%  achieve a perfect separation between the classes).
%
%  Area under ROC Curve (AUC) measure is very similar to Wilcoxon Rank Sum
%  Test and Mann-Whitney U Test. If "alg" parameter is set to 'Wilcoxon'
%  than AUC is calculated by performing that test and than canverting
%  result
%
%  The main properties of this code are:
%   * Ability to work with multi-dimensional data (X can have many columns).
%   * Ability to work with multi-class datasets ("y" can have more
%     than 2 different values).
%   * Speed - this code was written to calculate AUC's of large number of
%     features, fast.
%   * Returned AUC is always bigger than 0.5, which is equivalent of
%     testing for each feature colAUC(x,y) and colAUC(-x,y) and
%     returning the value of the bigger one.
%
% OUTPUT:
%   For 2 class problem colAUC returms a vector of AUC values one for each
%   feature/column.
%
%   For multi class problem AUC will be calculated for all combinations of
%   labels. All n!/((n-2)! 2!) = n(n-1)/2 of them (see nchoosek(n,2)),
%   where n is number of unique labels in "y" array.  Each pairing will be
%   returned in a separate row or "auc" and output "labels" will identify
%   the pair of labels.
%
%   For multi-class AUC "Total AUC" as defined by Hand & Till (2001) can be
%   calculated by mean(auc)}.
%
% See Also
%  * http://en.wikipedia.org/wiki/Receiver_operating_characteristic
%  * colAUC function from caTools package in R language
%    http://cran.r-project.org/web/packages/caTools
%
% Written by Jarek Tuszynski, SAIC, jaroslaw.w.tuszynski_at_saic.com
% Code covered by BSD License

%% define defaults and read input parameters
abs_auc = true;
alg  = 'Wilcoxon';
Plot = false;
nVarargs = length(varargin);
if mod(nVarargs,2)==1
  alg  = varargin{1};
end
for k = 1:nVarargs
  switch lower(varargin{k})
    case 'algorithm'
      alg = varargin{k+1};
    case 'plot'
      Plot = varargin{k+1};
    case 'abs'
      abs_auc = varargin{k+1};
  end
end

%% make sure inputs are in correct format
if (nargout==0), Plot = true; alg = 'ROC'; end
switch alg
  case 'ROC'
    alg = 'ROC';
  case {'Wilcoxon', 'Mann–Whitney U', 'Mann–Whitney–Wilcoxon', 'ranksum'}
    alg = 'Wilcoxon';
  otherwise
    alg = 'Wilcoxon';
end

%% Prepare for calculations & error-check
y = y(:);
[nr, nc] = size(X);                  % get dimentions of the data set
[uL, ~, yy] = unique(y);             % find all the classes among the labels
uL   = uL';
nL   = length(uL);                   % number of unique classes
if (nL<=1)
  error('colAUC: List of labels ''y'' have to contain at least 2 class labels.')
end
if (~isnumeric(X)), error('colAUC: ''X'' must be numeric'); end
if (nr~=length(y)), error('colAUC: length(y) and size(X,1) must be the same'); end
per  = nchoosek(1:nL,2);             % find all possible pairs of L columns
labels = uL(per);                    % those are pairs of labels stored as strings
if iscell(y)                         % unify label formats
  labelStr = uL;                     % string version of the labels
  y  = yy;
  uL = 1:nL;                         % numeric version of the labels
else
  labelStr = cell(size(uL));         % string version of the labels
  for i=1:nL, labelStr{i} = sprintf('%i',uL(i)); end
end
labelStr  = labelStr(per);           % those are pairs of labels stored as strings
L        = uL(ones(nr,1),:);
np       = size(per,1);              % how many possible pairs were found?
auc      = zeros(np,nc)+0.5;         % Initialize array to store results
rowLabel = cell(1,np);

%% prepare the plot, if needed
if (Plot)
  clf;
  hold on;
  color = 'bgrmcyk';
end

switch alg
  case 'ROC'
    %% Calculate AUC by integrating ROC curves
    iPlot = 0;
    for c=1:nc                           % for each column representing a feature
      [b, IDX] = sort( X(:,c));          % sort all columns and store them in X. IDX stores original positions
      nunq = find(diff(b)==0);           % find non-unique X's in column c (if vector is [1 1] nunq=1
      nTies = length(nunq);              % number of non-unique values
      if (nTies<nr-1)                    % make sure all numbers in X column are not the same
        IDX = y(IDX);                    % reorder label vector in the same order as X, or associate label with each number in X
        % assign column for each label (class) and for each point add 1 in the column corresponding to its class
        d = ( IDX(:,ones(1,nL)) == L );
        d = cumsum(d,1);                 % cumulative sum of columns or left node counts per class for all possible thresholds
        if(nTies), d(nunq, :) = []; end  % remove non unique rows (using nunq) if any
        d = [ zeros(1, nL); d ];         %#ok<AGROW> % append row of zeros at the beggining
        % assume that col#1 ploted on x axis is correct clasification and col#2 (y) is false find AUC
        for i=1:np                       % go through all permutations of columns in d
          c1 = per(i,1);                 % and identify 2 classes to be compared
          c2 = per(i,2);
          n  = d(end,c1)*d(end,c2);      % normalize area to 1 at the maximum
          if (n>0)
            xx=d(:,c1)/d(end,c1);
            yy=d(:,c2)/d(end,c2);
            auc(i,c) = trapz(xx, yy);    % Trapezoidal numerical integration
          else
            xx=[0,1];
            yy=[0,1];
          end
          if (Plot)
            %if abs_auc
            %  if (2*auc(i,c)<1), xx=1-xx; yy=1-yy; end % if auc<0.5 than mirror it to the other side of 0.5
            %end
            cc = mod(iPlot, length(color))+1;
            iPlot = iPlot+1;
            plot(xx, yy, color(cc));
            rowLabel{i} = sprintf('%s vs. %s', labelStr{i,1}, labelStr{i,2});
          end
        end
      end
    end
    
  case 'Wilcoxon'
    idxL = cell(1,nL);
    for i = 1:nL, idxL{i} = find(y==uL(i)); end
    for c = 1:nc                    % for each column representing a feature
      for i = 1:np                  % go through all permutations of columns in d
        c1 = per(i,1);              % and identify 2 classes to be compared
        c2 = per(i,2);
        x1 = X(idxL{c1},c);
        x2 = X(idxL{c2},c);
        n1 = length(x1);
        n2 = length(x2);
        if (n1>0 && n2>0)
          r = avrRank([x1(:); x2(:)]);
          U = sum(r(1:n1)) - n1*(n1+1)/2; % Wilcoxon rank-sum test or Mann–Whitney U test
          auc(i,c) = U / (n1*n2);
        end
      end % end of 'for i' looplength(
    end % end of 'for j' loop
end
%if abs_auc
%  auc = 0.5 + abs(0.5-auc); % if any auc<0.5 than mirror it to the other side of 0.5 auc is a matrix
%end

%% finalize the plot, if needed
if (Plot)
  %plot([0,1], [0,1], 'k')
  xlabel({'probability of false alarm','(1-Specificity)'});
  ylabel({'probability of detection'  ,'(Sensitivity)'});
  xlim([0,1.01]);
  ylim([0,1.01]);
  title('ROC Curves');
  grid on;
  hold off;
  if (nc==1 && np<20)  % if too many curves than skip the labels
    legend(rowLabel, 'Location', 'SouthEast');
  end
end
