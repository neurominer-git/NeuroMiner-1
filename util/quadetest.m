function STATS = quadetest(x, xnames, ynames, filename, varargin)
%QUADETEST: Quade test for non parametric two way ANalysis Of VAriance.
%This function performs the Quade test to analyze unreplicated complete block
%designs.
%Dana Quade in 1979 proposed a test that is often more powerful than the
%Friedman test. It also eliminates block differences but weights the raw data
%indicate possibly more marked treatment effects. Whereas the Friedman test is
%basically an extension of the sign test, the Quade test is effectively an
%extension of the Wilcoxon signed rank test and is equivalent to it when the
%treatments are two.
%
% Syntax: 	STATS=quadetest(X,alpha)
%      
%     Inputs:
%           X - data matrix
%           ALPHA - significance level (default = 0.05).
%     Outputs:
%           Quade Statistic
%           Multiple comparisons (eventually)
%
%      Example: 
%
% x=[115 142 36 91 28; 28 31 7 21 6; 220 311 108 51 117; 82 56 24 46 33;...
% 256 298 124 46 84; 294 322 176 54 86; 98 87 55 84 25];
%
%           Calling on Matlab the function: quadetest(x)
%
%           Answer is:
%
% QUADE TEST FOR IDENTICAL TREATMENT EFFECTS:
% TWO-WAY BALANCED, COMPLETE BLOCK DESIGNS
% --------------------------------------------------------------------------------
% Number of observation: 35
% Number of blocks: 7
% Number of treatments: 5
% --------------------------------------------------------------------------------
% F-statistic approximation
% Quade test statistic W: 10.3788
% F=W df-num=4 df-denom=24 - p-value (2 tailed): 0.0001
% --------------------------------------------------------------------------------
%  
% POST-HOC MULTIPLE COMPARISONS
% --------------------------------------------------------------------------------
% Critical value: 35.6981
% Absolute difference among mean ranks
%      0     0     0     0     0
%     18     0     0     0     0
%     51    69     0     0     0
%     69    87    18     0     0
%     63    81    12     6     0
% 
% Absolute difference > Critical Value
%      0     0     0     0     0
%      0     0     0     0     0
%      1     1     0     0     0
%      1     1     0     0     0
%      1     1     0     0     0
%
%           Created by Giuseppe Cardillo
%           giuseppe.    cardillo-edta@poste.it
%
%           modified by Nikolaos Koutsouleris 12/2018
% To cite this file, this would be an appropriate format:
% Cardillo G. (2009). QUADETEST: Quade test for non parametric two way ANalysis Of VAriance
% http://www.mathworks.com/matlabcentral/fileexchange/25926

%Input Error handling
p = inputParser;
addRequired(p, 'x',@(x) validateattributes(x,{'numeric'},{'2d','real','finite','nonnan','nonempty'}));
addOptional(p, 'alpha',0.05, @(x) validateattributes(x,{'numeric'},{'scalar','real','finite','nonnan','>',0,'<',1}));
addOptional(p, 'verbose', false, @(x) validateattributes(x,{'numeric'},{'binary'}));

if ~exist('xnames','var') || isempty(xnames), 
    xnames = cellstr([repmat('Partition_',size(x,1),1) num2str((1:size(x,1))')])';
else
    xnames = regexprep(xnames,'-','_');
end

if ~exist('ynames','var') || isempty(ynames), 
    ynames = cellstr([repmat('Model_',size(x,2),1) num2str((1:size(x,2))')])';
else
    ynames = regexprep(ynames,'-','_');
end
parse(p,x,varargin{:});
%assert(all(x(:,2) == fix(x(:,2))),'Warning: all elements of column 2 of input matrix must be whole numbers')
alpha = p.Results.alpha;
verbose = p.Results.verbose;

[r,c]=size(x); %dimension of the input matrix
R=zeros(r,c); %preallocation
%For each block, compute the ranks
for I=1:r
    R(I,:)=tiedrank(x(I,:));
end
%Compute the range of each block and then rank them.
Q=tiedrank(range(x,2));
%Compute a modified version of the Friedman matrix
rij=(R-(c+1)/2).*repmat(Q,1,c);
Ti=sum(rij);
T2=sum(Ti.^2);
rij2=sum(sum(rij.^2));
T3=T2/r;
T4=rij2-T3;
k=r-1;
W=k*T3/T4; %The Quade statistic.
%The Quade statistic is approximable with the F distribution.
dfn=c-1;
dfd=dfn*k;
pvalue=1-fcdf(W,dfn,dfd);
t_blocks = table(r*c,r,c,'VariableNames',{'Observations','Blocks','Treatments'});
t_quade = table(W,dfn,dfd,pvalue,'VariableNames',{'W','DF_numerator','DF_denominator','two_tailed_p_value'});
%display results
if verbose
    tr=repmat('-',1,80); %set the divisor
    disp('QUADE TEST FOR IDENTICAL TREATMENT EFFECTS: TWO-WAY BALANCED, COMPLETE BLOCK DESIGNS')
    disp(tr)
    disp(t_blocks)
    disp('QUADE''S STATISTICS: F-statistic approximation')
    disp(tr)
    disp(t_quade)

end

if pvalue<alpha
   
    tmp=repmat(Ti,c,1); Rdiff=abs(tmp-tmp'); %Generate a matrix with the absolute differences among ranks
    % Critical value based on Heckert and Filliben (2003) which uses the student-t distribution.
    denom   = realsqrt(2*r*T4/dfd);
    cv      = tinv(1-alpha/2,dfd)*denom; %critical value
    mc      = Rdiff>cv; %Find differences greater than critical value
    tval    = Rdiff ./ denom;
    pval    = 1-tcdf(tval,dfd);
    % We do need to correct only the p values in the lower triangle of the
    % matrix without the diagnonal for multiple comparisons
    Jx       = itril(size(pval),0);
    Jy       = itriu(size(pval),1);
    pval_nan = pval; pval_nan(Jx) = nan;
    [mc_fdr, cv_p, ~,pval_fdr] = fdr_bh(pval_nan(~isnan(pval_nan)),0.05,'pdep');
    pval_nan_fdr = pval_nan;
    pval_nan_fdr(Jy) = pval_fdr;
    pval(Jx) = nan;
    mc_nan_fdr = nan(size(pval));
    mc_nan_fdr(Jy) = mc_fdr;
    %display results
    if verbose
        disp(' ')
        disp('POST-HOC MULTIPLE COMPARISONS')
        disp(tr)
        fprintf('Critical value: %0.4f\n',cv)
        disp('Absolute difference among mean ranks')
        %disp(tril(Rdiff))
        disp(tril(Rdiff))
        disp('Absolute difference > Critical Value')
        disp(tril(mc))
    end
end

if nargout >0
    STATS.alpha = alpha;
    STATS.T3 = T3;
    STATS.T4 = T4;
    STATS.W = W; 
    STATS.p = pvalue;
    STATS.tbl_blk = t_blocks;
    STATS.tbl_quade = t_quade;
    if exist('tmp','var')
        STATS.Rank_diffs = Rdiff;
        STATS.Crit_val = cv;
        STATS.Crit_met = mc;
        STATS.Tval = tval;
        STATS.Crit_pval_fdr = cv_p;
        STATS.Crit_met_fdr = mc_fdr;
        STATS.Pval = pval;
        STATS.Pval_fdr = pval_fdr;
        STATS.tbl_perf = array2table(x,'VariableNames',ynames', 'RowNames',xnames);
        STATS.tbl_rankdiff_posthoc = array2table(Rdiff,'VariableNames',ynames', 'RowNames',ynames);
        STATS.tbl_p_posthoc = array2table(pval,'VariableNames',ynames', 'RowNames',ynames);
        STATS.tbl_p_fdr_posthoc = array2table(pval_nan_fdr,'VariableNames',ynames', 'RowNames',ynames);
        STATS.tbl_crit_posthoc = array2table(mc,'VariableNames',ynames, 'RowNames',ynames');
        STATS.tbl_crit_fdr_posthoc = array2table(mc_nan_fdr,'VariableNames',ynames', 'RowNames',ynames);
        if exist('filename','var') && ~isempty(filename)
            sheetnames = { 'pp', 'blk', 'Quade', 'Rankdiff', 'Crit', 'P_uncorr','P_fdr' };
            tblnames = { 'tbl_perf', 'tbl_blk', 'tbl_quade', 'tbl_rankdiff_posthoc','tbl_crit_posthoc','tbl_p_posthoc', 'tbl_p_fdr_posthoc'} ;
            display_writetable(STATS,filename, tblnames, sheetnames)
        end
    end
else
    STATS=[];
end
