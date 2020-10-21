% =========================================================================
% Coded by Youngmin Ha (y.ha.1@research.gla.ac.uk) on 21 Jan 2015
% =========================================================================
% This program demonstrates multivariate relevance vector regression (MRVR).
%   A user can choose either origianl MRVR algorithm devloped by 
%   Thayananthan or fast MRVR algorithm devloped by Ha
% =========================================================================

clear, clc, close all;

%% free parameters
rng(7)      % seed of random number generation
isFast      = true; % true (fast MRVR) or false (origianl MRVR)
N           = 200;  % # of training smaples
maxIts      = 1000; % maximum iteration number of EM algorithm
tolerance   = .1;   % tolerance value to check convergence of EM algorithm
distX       = 'uniform'; % distribution of training samples in x-axis

% kernel
kernelType  = '+gauss';
kernelWidth	= 1.6;

%% generate training samples
% x
switch distX
    case 'uniform'
        x = 20.*(rand(N,1)-.5);
    case 'gauss'
        x = 3*randn(N,1);
    otherwise
        error('unsupported distX')
end
x = sort(x);
xLin = linspace(-10,10,N)';

% generate sinc function
y1 = sinc((x)./pi);
y1lin = sinc((xLin)./pi);

% generate linear function
y2 = 0.05.*x;
y2lin = 0.05.*xLin;

% multivariate output
Y =[y1 y2];
Ylin = [y1lin y2lin];

% add Gaussian random noise
noise = rand(1,2); % standard deviation of noise
rho = -1 + 2.*rand; % correlation coefficient of noise
R = eye(2); % R is correlation matrix
R(1,2) = rho;
R(2,1) = rho;
Omega = diag(noise)*R*diag(noise); % 2 by 2 random covariance matrix
E = mvnrnd([0 0],Omega,N);
T = Y + E;

%% multivariate relevance vector regression
Phi	= sbl_kernelFunction(x,x,kernelType,kernelWidth);

tic
if isFast
    [used, ~, Mu, invSigma, OmegaHat] = fmrvr(Phi,T,maxIts,tolerance);
else
    [used, ~, Mu, invSigma, OmegaHat] = mrvr(Phi,T,maxIts,tolerance);
end
toc

Phi	= sbl_kernelFunction(xLin,x,kernelType,kernelWidth);
Ypredict = Phi(:,used)*Mu;
varNoise = repmat(diag(OmegaHat)',N,1);
if isFast
    stdPredict = sqrt(varNoise.*repmat(1+diag(Phi(:,used)...
        *(invSigma\Phi(:,used)')),1,2));
else
    for j = 1:2
        varWeight(:,j) = diag(Phi(:,used)*(invSigma(:,:,j)\Phi(:,used)'));
    end
    stdPredict = sqrt(varNoise + varWeight);
end

% check if bias is used or not
if kernelType(1) == '+'
    used = used - 1;
    if used(1) == 0
        used(1)	= [];
    end
end

%% plot
fig = figure('Position',[100 50 500 900]); % [left, bottom, width, height]
set(fig,'defaulttextinterpreter','latex');
for j = 1:2
    subplot(2,1,j)
    hold on;
    plot(x,T(:,j),'black.','MarkerSize',14)
    h_y = plot(xLin,Ylin(:,j),'r:','LineWidth',2);
    h_yPredict = plot(xLin,Ypredict(:,j),'b-','LineWidth',2);
    h_pm = plot(xLin,Ypredict(:,j) + stdPredict(:,j),'b--','LineWidth',2);
    plot(xLin,Ypredict(:,j) - stdPredict(:,j),'b--','LineWidth',2);
    h_rv = plot(x(used),T(used,j),'go','LineWidth',2,'MarkerSize',8);
    
    subJ = num2str(j);
    switch j
        case 1
            legend([h_yPredict h_pm],...
                {['Predicted mean ' '$(y_{*,' subJ '})$'],...
                ['$y_{*,' subJ '}\pm\sigma_{*,' subJ '}$']},...
                'Location','Best','Interpreter','Latex')
        case 2
            legend([h_y h_yPredict h_pm h_rv],...
                {'True function',...
                ['Predicted mean ' '($y_{*,' subJ '}$)'],...
                ['$y_{*,' subJ '}\pm\sigma_{*,' subJ '}$'],...
                'Relevance vectors'},...
                'Location','Best','Interpreter','Latex')
    end
    xlabel('$x$')
    ylabel({['$y_' subJ '$']})
    if j == 1
        title('1st output')
    else
        title('2nd output')
    end
    hold off
end

%% RMS and estimated noise level
for j = 1:2
    switch j
        case 1
            text = 'st';
        case 2
            text = 'nd';
    end
    fprintf(['%d' text ' output\n'],j)
    fprintf('Regression error (RMS): %g\n', ...
	sqrt(mean((Y(:,j)-Ypredict(:,j)).^2)))
    fprintf('Estimated noise level: %.4f (true: %.4f)\n\n', ...
        sqrt(OmegaHat(j,j)), noise(j))
end

%% measure the performance of covariance matrix estimation
V = size(Y,2); % V is output dimensionality
temp = (Omega\OmegaHat)'; % Omega is true, and OmegaHat is estimated
loss1 = trace(temp) - log(det(temp)) - V; % entropy loss
loss2 = trace((temp - eye(V))^2); % quadratic loss

disp(['Entropy loss = ' num2str(loss1)])
disp(['Quadratic loss = ' num2str(loss2)])