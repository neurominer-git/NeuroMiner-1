
% -------------------------------------------------------------------------
% Function: out=preprocessing(X,method)
% -------------------------------------------------------------------------
% Aim:
% Data preprocessing
% ------------------------------------------------------------------------
% Input: 
% X, matrix (n,p), data matrix n-objects, p-variables 
% method, string defining type of data preprocessing
% 
% classical preprocessing:
% 'mean centering' - columnwise centering of X
% 'standardisation' - columnwise standardisation of X
% 'SNV' - Standard Normal Variate
%
% robust preprocessing:
% 'median centering' - columnwise centering of X with classical median
% 'l1-median centering' - columnwise centering of X with L1-median
% 'qn-standardisation' - columnwise standardisation of X with Qn-estimator
% 'qn-autoscaling' - autoscaling by centering around L1-median and Qn
% standardisation
% 'sn-standardisation' - columnwise standardisation of X with Sn-estimator
% 'sn-autoscaling' - autoscaling by centering around L1-median and Sn
% standardisation
% 'mad' - columnwise Median of Absolute Deviation
% 'median-snv' - Standard Normal Variate with classical median
% 'sn-snv' - Standard Normal Variate with Qn-estimator
% 'sn-snv' - Standard Normal Variate with Sn-estimator
% ------------------------------------------------------------------------
% Output:
% out, preprocessed X
% -----------------------------------------------------------------------
% Example: 
% out=preprocessing(X,'mean centering')

% Written by Sven Serneels
% MiTAC, University of Antwerp
% December 2004

function out = preprocessing(X,method)

% Check if 'method' string is peoperly specified
options=[{'mean centering'};{'standardisation'};{'snv'};{'median centering'};{'l1-median centering'};,...
        {'qn-standardisation'};{'qn-autoscaling'};{'sn-standardisation'};{'sn-autoscaling'};{'mad'};,...
        {'median-snv'};{'qn-snv'};{'sn-snv'}];

if isempty(strmatch(lower(method),lower(options),'exact'))
    errordlg([{'Input variable ''method'' is wrongly specified'};{' '};{'See: help preprocessing'}] ,'Error')
    out=[]
    return
end

[m,n]=size(X);

switch method
    
    % Classical proprocessing
    case 'mean centering'
        out=X-(ones(m,1)*mean(X));
        
    case 'standardisation'
        out=X./(ones(m,1)*std(X));
        
    case 'snv'
        out=(X-mean(X')'*ones(1,n))./(std(X')'*ones(1,n));
        
    % Robust proprocessing
    case 'median centering'
        out=X-(ones(m,1)*median(X));
        
    case 'l1-median centering'
        out=X-(ones(m,1)*L1median(X));
        
    case 'qn-standardisation'
        out=X./(ones(m,1)*qn(X));
        
    case 'qn-autoscaling'
        out=X-(ones(m,1)*L1median(X));
        out=out./(ones(m,1)*qn(out));
        
    case 'sn-standardisation'
        out=X./(ones(m,1)*sn(X));
              
    case 'sn-autoscaling'
        out=X-(ones(m,1)*L1median(X));
        out=out./(ones(m,1)*sn(out));
        
    case 'mad'
        out=mad(X);
        
    case 'median-snv'
        out=(X-median(X')'*ones(1,n))./(median(X')'*ones(1,n));
        
    case 'qn-snv'
        out=(X-qn(X')'*ones(1,n))./(qn(X')'*ones(1,n));
        
    case 'sn-snv'
        out=(X-sn(X')'*ones(1,n))./(sn(X')'*ones(1,n));
end
    

% ---> Median Absolute Deviation
function m=mad(X)

% MAD computes the median absolute deviation of X. If X is
%  a matrix, MAD is a row vector containing the MAD's of the
%  columns of X.
%
% ! Includes correction for consistency !
%
% Written by S. Serneels, 17.12.2003

[n,p]=size(X);
Xmc=X-repmat(median(X),n,1);
m=1.482*median(abs(Xmc));


% ---> QNS
function s=qn(y)

% QN Qn scale estimator
% ------------------------------------------
% Input:  y, matrix of size (n,p)
% ------------------------------------------
% Output: s, vector of size (1,p) containing the Qn scale estimates of the
% columns of y
% ------------------------------------------
% The Qn estimator is proposed in P.J. Rousseeuw, C. Croux, Alternatives to
% the median absolute deviation, J. Am. Statist. Assoc., 88 (1993),
% 1273-1283

% Written by Sven Serneels, University of Antwerp

if size(y,2)>1
    if size(y,1)>1
        for i=1:size(y,2)
            s(:,i)=qnsven(y(:,i));
        end
    else
        y=y';
        s=qnsven(y);
    end
else
    s=qnsven(y);
end;


% ---> Sn scale estimator
function s=sn(y)

% SN Sn scale estimator
% ------------------------------------------
% Input:  y, matrix of size (n,p)
% ------------------------------------------
% Output: s, vector of size (1,p) containing the Sn scale estimates of the
% columns of y
% ------------------------------------------
% The Sn estimator is proposed in P.J. Rousseeuw, C. Croux, Alternatives to
% the median absolute deviation, J. Am. Statist. Assoc., 88 (1993),
% 1273-1283

% Written by Sven Serneels, University of Antwerp

if size(y,2)>1
    if size(y,1)>1
        for i=1:size(y,2)
            s(:,i)=snsven(y(:,i));
        end
    else
        y=y';
        s=snsven(y);
    end
else
    s=snsven(y);
end;

% -----------------------------------------
function s=snsven(y)

n=length(y);
if n>1000
    sy=sort(y);
    nbins=floor(n/10);
    mys=zeros(nbins,1);
    ninbins=floor(n/nbins);
    for i=1:nbins
        if (mod(n,nbins)~=0 & i==nbins)
            mys(i)=median(sy((i-1)*ninbins+1:n));
        else
            mys(i)=median(sy((i-1)*ninbins+1:i*ninbins));
        end
    end
    y=mys;
    n=nbins;
end
pairwisediff=sort(abs(repmat(y',n,1)-repmat(y,1,n))); 
pairwisediff=pairwisediff(floor((n+1)/2),:);
pairwisediff=sort(pairwisediff);
s=1.1926*(pairwisediff(floor(n/2)+1));

cn=1;
switch n
    case 2 
        cn=0.743;
    case 3 
        cn=1.851;
    case 4 
        cn=0.954;
    case 5
        cn=1.351;
    case 6 
        cn=0.993;
    case 7 
        cn=1.198;
    case 8 
        cn=1.005;
    case 9 
        cn=1.131;
    otherwise 
        if (mod(n,2)==1) 
            cn=n/(n-0.9);
        end
end
s=cn*s;


%-----------------------------
function [mX]=L1median(X,tol);

% L1MEDIAN calculates the multivariate L1-median 
% I/O: [mX]=L1median(X,tol);
%
% X is the data matrix 
% tol is the convergence criterium; the iterative proces stops when ||m_k - m_{k+1}|| < tol.
%
% Ref: Hossjer and Croux (1995) "Generalizing Univariate Signed Rank Statistics for Testing
% and Estimating a Multivariate Location Parameter", Non-parametric
% Statistics, 4, 293-308

if nargin <2
    tol=1.e-08;
end;

[n,p]=size(X);
maxstep=100;

% initializing starting value for m
m=median(X);
k=1;
while (k<=maxstep)
    mold=m;
    Xext=sortrows([norme(X-repmat(m,n,1)) X],1);
    dx=Xext(:,1);
    X=Xext(:,2:p+1);
    if all(dx)
        w=1./dx;
    else
        ww=dx(all(dx,2));
        w=1./ww;
        w=[zeros(length(dx)-length(w),1);w];
    end
    delta=sum((X-repmat(m,n,1)).*repmat(w,1,p),1)./sum(w);
    nd=norme(delta);
    if all(nd<tol)
        maxhalf=0;
    else
        maxhalf=log2(nd/tol);
    end
    m=mold+delta;   % computation of a new estimate
    nstep=0;
    while all(mrobj(X,m)>=mrobj(X,mold))&(nstep<=maxhalf)
        nstep=nstep+1;
        m=mold+delta./(2^nstep);
    end
    if (nstep>maxhalf)
        mX=mold;
        break
    end
    k=k+1;
end

mX=m;

%-----
function n=norme(X);

% NORME calculates the euclidian norm of matrix X
% the output is a columnvector containing the norm of each row
% I/O: n=norme(X);

n=sqrt(sum(X.^2,2));

%--------
function s=mrobj(X,m)

% MROBJ computes objective function in m based on X and a

xm=norme(X-repmat(m,size(X,1),1));
s=sum(xm,1)';


%-------------------------------------------
function s=qnsven(y)

n=length(y);
% Do binning for big n
if n>1000
    sy=sort(y);
    nbins=floor(n/10);
    mys=zeros(nbins,1);
    ninbins=floor(n/nbins);
    for i=1:nbins
        if (mod(n,nbins)~=0 & i==nbins)
            mys(i)=median(sy((i-1)*ninbins+1:n));
        else
            mys(i)=median(sy((i-1)*ninbins+1:i*ninbins));
        end
    end
    y=mys;
    n=nbins;
end
h=floor(n/2)+1;
k=0.5*h*(h-1);
pairwisediff=repmat(y,1,n)-repmat(y',n,1); 
pairwisediff=sort(abs(pairwisediff(find(tril(ones(n,n),-1)))));
s=2.2219*(pairwisediff(k));


switch n
    case 1
        error('Sample size too small');
    case 2
        dn=0.399;
    case 3 
        dn=0.994;
    case 4 
        dn=0.512;
    case 5 
        dn=0.844;
    case 6 
        dn=0.611;
    case 7 
        dn=0.857;
    case 8 
        dn=0.669;
    case 9 
        dn=0.872;
        
    otherwise 
        if (mod(n,2)==1) 
            dn=n/(n+1.4);
        elseif (mod(n,2)==0) 
            dn=n/(n+3.8);
        end
end
s=dn*s;