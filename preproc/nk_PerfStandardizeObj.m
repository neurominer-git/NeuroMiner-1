function [sY, IN] = nk_PerfStandardizeObj(Y, IN)
% =========================================================================
% FORMAT function [sY, IN] = nk_PerfStandardizeObj(Y, IN)
% =========================================================================
% Standardisation of matrix Y using different types of estimators of mean 
% and spread of each feature (in sIND). Standardization can be applied to 
% a given destination subgroup (logical index dIND) and estimators can be 
% computed from a given source soubgroup (logical index sIND).
%
% For standardization using median or mean:
% If WinsOpt is provided than any values outside the WinsOpt are clamped to
% the WinsOpt limits. In case of Winsorization the standardized data is 
% re-centered so that the mean of each column is again zero (in this case
% meanY2 is returned to the parent function).
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Nikolaos Koutsouleris, 04/2018

% =========================== WRAPPER FUNCTION ============================ 
if iscell(Y) && exist('IN','var') && ~isempty(IN)
    sY = cell(1,numel(Y)); 
    for i=1:numel(Y), 
        % Define active indices depending on training or testing situation
        if isfield(IN,'trained') && IN.trained 
           if isfield(IN,'sTsInd'), IN.sIND = IN.sTsInd{i}; else IN.sIND =[]; end
           if isfield(IN,'dTsInd'), IN.dIND = IN.dTsInd{i}; else IN.dIND = []; end
        else
           if isfield(IN,'sTrInd'), IN.sIND = IN.sTrInd{i}; else IN.sIND = []; end
           if isfield(IN,'dTrInd'), IN.dIND = IN.dTrInd{i}; else IN.dIND =[]; end
        end
        sY{i} = PerfStandardizeObj(Y{i}, IN ); 
    end
else
    if ~exist('IN','var'), IN=[]; end
    % Define active indices depending on training or testing situation
    if isfield(IN,'trained') && IN.trained
       if isfield(IN,'sTsInd'), IN.sIND = IN.sTsInd; else IN.sIND =[]; end
       if isfield(IN,'dTsInd'), IN.dIND = IN.dTsInd; else IN.dIND =[]; end
    else
       if isfield(IN,'sTrInd'), IN.sIND = IN.sTrInd; else IN.sIND = []; end
       if isfield(IN,'dTrInd'), IN.dIND = IN.dTrInd; else IN.dIND = []; end
    end
    [ sY, IN ] = PerfStandardizeObj(Y, IN );
end
% ========================================================================= 
function [sY, IN] = PerfStandardizeObj(Y, IN)

global VERBOSE
[mY, nY] = size(Y);

if ~isfield(IN,'method'),   IN.method = 'standardization using median'; end
if ~isfield(IN,'trained'),  IN.trained = false; end

if isempty(IN), eIN=true; else, eIN=false; end
% Zero-out features with non-finite values
if eIN || ~isfield(IN,'zerooutflag') || isempty(IN.zerooutflag), IN.zerooutflag = 2;  end
if eIN || ~isfield(IN,'WINSOPT') || isempty(IN.WINSOPT), IN.WINSOPT= [];  end

% Compute standardization parameters from training data if not available in
% IN ==> mean (std) for each feature
if eIN || ~IN.trained
    
    % with sIND the user can specify which subgroup of cases needs to be
    % used for computing the Z-transformation params
    if ~isfield(IN,'sIND') || isempty(IN.sIND), IN.sIND = true(mY,1); end
    if ~islogical(IN.sIND), IN.sIND = logical(IN.sIND); end 
    nG = size(IN.sIND,2);
    
    switch IN.method
        
     
        case {'standardization using median', 'standardization using mean'}
            
            funmean = 'nm_nanmedian'; if strcmp(IN.method,'standardization using mean'), funmean = 'nm_nanmean'; end
            IN.meanY = zeros(nG,nY);
            IN.stdY = zeros(nG,nY);

            for i=1:nG    
                if ~sum(IN.sIND(:,i)), continue, end
                IN.meanY(i,:) = feval(funmean,(Y(IN.sIND(:,i),:))); % compute mean
                nanvec = isnan(IN.meanY(i,:));
                if any(nanvec), fnanvec = find(nanvec); IN.meanY(i,fnanvec) = feval(funmean,Y(IN.sIND(:,i),fnanvec)); end
                IN.stdY(i,:) = nm_nanstd(Y(IN.sIND(:,i),:)); % and STD
                std_zeros = find(IN.stdY(i,:)<=realmin);
                if ~isempty(std_zeros)
                    cprintf('red','\tWarning: %g non-discriminative feature(s) exist in the data! ', numel(std_zeros)); cprintf('black',' ');
                    IN.stdY(i,std_zeros) = 10e-6;
                end
            end

            if isfield(IN,'WINSOPT') && ~isempty(IN.WINSOPT) && IN.WINSOPT
                IN.meanY2 = zeros(nG,nY);
            end
        % X_i = (X_i - mean(X)) 
        case 'mean-centering'
            funmean = 'nm_nanmean';
            IN.meanY = zeros(nG,nY);
            for i=1:nG    
                if ~sum(IN.sIND(:,i)), continue, end
                IN.meanY(i,:) = feval(funmean,(Y(IN.sIND(:,i),:))); % compute mean
                nanvec = isnan(IN.meanY(i,:));
                if any(nanvec), fnanvec = find(nanvec); IN.meanY(i,fnanvec) = nm_nanmean(Y(IN.sIND(:,i),fnanvec)); end
            end
        case 'l1-median centering' 
            IN.L1median = zeros(nG,nY); for i=1:nG, IN.L1median(i,:) = L1median(Y(IN.sIND(:,i),:)); end
        case 'qn-standardization'
            IN.qn = zeros(nG,nY); for i=1:nG, IN.qn(i,:) = qn(Y(IN.sIND(:,i),:)); end
        case 'sn-standardization'
            IN.sn = zeros(nG,nY); for i=1:nG, IN.sn(i,:) = sn(Y(IN.sIND(:,i),:)); end
    end
    IN.trained = 1;
    
else
    switch IN.method
        case {'standardization using median', 'standardization using mean', 'mean-centering'}
            nG = size(IN.meanY,1);
        case 'l1-median centering'
            nG = size(IN.L1median,1);
        case 'qn-standardization'
            nG = size(IN.qn,1);
        case 'sn-standardization'
            nG = size(IN.sn,1);
    end
    if isempty(IN.sIND), IN.sIND = true(mY,1); end
end
% with dIND the user defined which cases in the matrix are normalized 
if eIN ||~isfield(IN,'dIND')  || isempty(IN.dIND), 
    if size(IN.sIND,2) > 1,
        error('Destination [%g-by-%g] matrix is missing in IN structure', mY, nG);
    else
        IN.dIND = true(mY,1);
    end
end
if ~islogical(IN.dIND), IN.dIND = logical(IN.dIND); end 

sY = zeros(size(Y)); 
% Standardize (sub)matrix

for i=1:nG
    
    if ~sum(IN.dIND(:,i)), continue, end
    % Perform standardization
    switch IN.method
        case {'standardization using median', 'standardization using mean'}
            funmean = 'nm_nanmedian'; if strcmp(IN.method,'standardization using mean'), funmean = 'nm_nanmean'; end
            tY = bsxfun(@rdivide, bsxfun(@minus, Y(IN.dIND(:,i),:), IN.meanY(i,:)), IN.stdY(i,:));
            % Winsorize data to prespecified highest and lowest values
            if ~isempty(IN.WINSOPT) && ~sum(IN.meanY2(:))
                %Perform standardization
                indP = tY>IN.WINSOPT; indN = tY<-1*IN.WINSOPT;
                tY(indP) = IN.WINSOPT; tY(indN) = -1*IN.WINSOPT;
                % Re-center the data to the mean
                IN.meanY2(i,:) = feval(funmean,(tY(IN.sIND(:,i),:)));
                nanvec = isnan(IN.meanY2(i,:));
                if any(nanvec), fnanvec = find(nanvec); IN.meanY2(i,fnanvec) = nm_nanmean(Y(IN.sIND(:,i),fnanvec));end
                tY = bsxfun(@minus, tY, IN.meanY2(i,:));
                if VERBOSE, fprintf('\n\t\t\tWinsorized # of values >%gSD threshold: (+) %g / (-) %g; data re-centered.', IN.WINSOPT, sum(indP(:)),sum(indN(:))); end
            end
        case 'mean-centering'
            tY = bsxfun(@minus, Y(IN.dIND(:,i),:), IN.meanY(i,:));
        case 'l1-median centering'
            tY = bsxfun(@minus, Y(IN.dIND(:,i),:), ones(size(Y(IN.dIND(:,i),:),1),1)*IN.L1median(i,:));
        case 'qn-standardization'
            tY = bsxfun(@rdivide, Y(IN.dIND(:,i),:), ones(size(Y(IN.dIND(:,i),:),1),1)*IN.qn(i,:));
        case 'sn-standardization'
            tY = bsxfun(@rdivide, Y(IN.dIND(:,i),:), ones(size(Y(IN.dIND(:,i),:),1),1)*IN.sn(i,:));
    end

    sY(IN.dIND(:,i),:) = tY;
    
end

[ sY, IN ] = nk_PerfZeroOut(sY, IN);
% Important: copy original matrix Y to sY!

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

function n=norme(X)

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
