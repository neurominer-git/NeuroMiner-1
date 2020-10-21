% =========================================================================
% Coded by Youngmin Ha (y.ha.1@research.gla.ac.uk) on 21 Jan 2015
%   by referring to Ha, Youngmin. "Fast multivariate relevance vector 
%   regression," submitted to Annals of Mathematics and Artificial 
%   Intelligence (2015).
% =========================================================================
% FUNCTION fmrvr performs fast multivariate relevance vector regression
%   developed by Ha.
% ******************************** Input **********************************
% Phi:      matrix, design matrix
% T:        column vector or matrix, target (i.e. true function + noise)
% maxIters:	scalar, maximum iteration number of EM algorithm
% tolerance:scalar, tolerance value to check convergence of EM algorithm
% ******************************** Output *********************************
% used:     column vector, index of relevance vectors
% alpha: 	column vector, inverse variance of weight
% Mu:       column vector or matrix, mean of weight
% invSigma: matrix, inverse of covariance matrix of weight
% Omega:    scalar or matrix, estimated covariance matrix of noise
% nIters:   scalar, number of iterations of EM algorithm
% =========================================================================

function [used, alpha, Mu, invSigma, Omega, nIters] = fmrvr(Phi, T, ...
    maxIters, tolerance)
%% constants
N = size(T,1); % # of training samples
V = size(T,2); % # of output dimensions

assert(size(Phi,1) == N, 'unexpected matrix size of Phi')
assert(size(Phi,2) == N + 1, 'unexpected matrix size of Phi')

%% initialisation
Omega = .1.*cov(T);
alpha = inf(N+1,1);
mask = logical(alpha < inf); % false means alpha=inf, true means alpha<inf

%% EM algorithm
isConverged = false;
for iterNum = 1:maxIters
    % update s', q', s, and q
    sPrime = nan(N+1,1);
    qPrime = nan(N+1,V);
    s = nan(N+1,1);
    q = nan(N+1,V);
    allAlphaInf = all(~mask);
    for i = 1:N+1
        if allAlphaInf % if all alpha are inf
            sPrime(i) = Phi(:,i)'*Phi(:,i);
            qPrime(i,:) = Phi(:,i)'*T;
        else
            sPrime(i) = Phi(:,i)'*Phi(:,i) - Phi(:,i)'*PhiSigmaPhi*Phi(:,i);
            qPrime(i,:) = Phi(:,i)'*T - Phi(:,i)'*PhiSigmaPhi*T;
        end
        
        if ~mask(i) % isinf(alpha(i))
            s(i) = sPrime(i);
            q(i,:) = qPrime(i,:);
        else
            temp = alpha(i)./(alpha(i) - sPrime(i));
            s(i) = temp.*sPrime(i);
            q(i,:) = temp.*qPrime(i,:);
        end        
    end
    
    % select i which maximizes deltaL(i)
    deltaL = nan(N+1,1);
    task = repmat('non', N+1, 1);
    theta = zeros(N+1,1);
    alphaNew = nan(N+1,1);
    invOmega = inv(Omega);
    for i = 1:N+1
        theta(i) = trace(invOmega*(q(i,:)'*q(i,:)))./V - s(i);
        qPrimeSq = trace(invOmega*(qPrime(i,:)'*qPrime(i,:)));
        if theta(i) > 0 % if i should be in model
            if mask(i) % ~isinf(alpha(i)) % if i was already in model
                % reestimate alpha_i
                task(i,:) = 'est';
                alphaNew(i) = (s(i).^2)./theta(i);
                % if alphaNew is 0, deltaL is -inf
                if alphaNew(i) ~= 0
                    temp = 1./alphaNew(i) - 1./alpha(i);
                    deltaL(i) = qPrimeSq./(sPrime(i) + 1./temp) ...
                        - V.*log(1 + sPrime(i).*temp);
                else
                    deltaL(i) = -inf;
                end
            else
                % add i
                task(i,:) = 'add';
                temp = qPrimeSq./sPrime(i);
                deltaL(i) = temp - V + V.*log(V./temp);
            end
        elseif mask(i) % ~isinf(alpha(i)) % if i was already in model
            % delete i
            task(i,:) = 'del';
            deltaL(i) = qPrimeSq./(sPrime(i) - alpha(i)) ...
                - V.*log(1 - sPrime(i)./alpha(i));
        end
        if imag(deltaL(i))
            % deltaL(i) has a imaginary number if inside of log is negative,
            %   but it is numerical error. Therefore, replace it with -inf
            deltaL(i) = -inf;
            warning('deltaL has an imaginary number')
        end
    end
    if all(isnan(deltaL)) % if all tasks are non
        break
    elseif all(deltaL(~isnan(deltaL)) == -inf)
        warning('all values of deltaL are -inf')
        i = randi(N+1);
    else
        [~, i] = max(deltaL);
    end
    
    % update alpha(i)
    switch task(i,:)
        case 'est'
            changeLogAlpha = log(alpha(i)./alphaNew(i));
            alpha(i) = alphaNew(i);
            if abs(changeLogAlpha) < tolerance
                if all(theta(~mask) <= 0) % ~mask means "out of the model"
                    isConverged = true;
                end
            end
        case 'add'
            alpha(i) = (s(i).^2)./theta(i);
            mask(i) = true;
        case 'del'
            alpha(i) = inf;
            mask(i) = false;
        otherwise
            error('task should be one of est, add, and del')
    end
    
    % update Omega
    if iterNum ~= 1
        Omega = T'*(T - PhiMask*Mu)./N;
    end

    % update PhiMask, invSigma, and Mu
    PhiMask = Phi(:,mask);
    invSigma = PhiMask'*PhiMask + diag(alpha(mask));
    assert(all(all(invSigma == invSigma')), 'invSigma should be symmetric')
    SigmaPhi = invSigma\PhiMask';
    Mu = SigmaPhi*T;
    PhiSigmaPhi = PhiMask*SigmaPhi;
        
    if isConverged
        break
    end
end
nIters = iterNum;
disp(['2. # of iterations = ' num2str(iterNum)])

% index of relevance vectors
used = find(mask);