% =========================================================================
% Coded by Youngmin Ha (y.ha.1@research.gla.ac.uk) on 21 Jan 2015
%   by referring to Ha, Youngmin. "Fast multivariate relevance vector 
%   regression," submitted to Annals of Mathematics and Artificial 
%   Intelligence (2015).
% =========================================================================
% This program conducts Monte Carlo simulations to compare between
%   origianl multivariate relevance vector regression (MRVR) algorithm 
%   devloped by Thayananthan and fast MRVR algorithm devloped by Ha
% =========================================================================

clear, clc, close all

%% free parameters
nSims   = 101; % # of simulations
maxIts  = 1000; % maximum iteration of EM algorithm
tolerance = .1; % tolerance value to check convergence of EM algorithm
N       = round(linspace(50,300,6)); % # of training samples
V       = 1:5; % # of output dimensions
shiftX  = [0 -1 1 -2 2].*2; % shift in x-axis of true functions

% kernel
kernelType = '+gauss';
kernelWidth = 1.6;

%% memory allocation
nIters  = nan(nSims,length(V),length(N),2);
runTime = nIters;
loss1   = nIters;
loss2   = nIters;
rmse    = nIters;
nRvs    = nIters;

%% loop of MC simulations
for i = 1:nSims
    rng(i) % seeds the random number generator for each simulation
    for j = 1:length(V)
        isCovPos = false;
        
        % generate a ranodm positive semi-definite matrix
        while isCovPos == false
            % R is correlation matrix of noise
            R = tril(-1 + 2.*rand(V(j)), -1);
            assert(all(all(R ~= -1 & R ~= 1)), 'covariance matrix is singular')
            R = R + R';
            R(logical(eye(V(j)))) = ones(V(j),1);

            % Omega is covariance matrix of noise
            Sigma = rand(1,V(j)); % standard deviation of noise
            Omega = diag(Sigma)*R*diag(Sigma);

            % check covariance matrix is positive semi-definite
            try
                [~] = mvnrnd(zeros(1,V(j)), Omega);
                isCovPos = true;
            catch
                if(strfind(lasterr, 'SIGMA must be a symmetric positive semi-definite matrix'))
                    isCovPos = false;
                end
            end
        end

        for k = 1:length(N); % different # of training samples
            x = linspace(-10, 10, N(k))';

            % multivariate output
            Y = [];
            for l = 1:V(j)
                Y = [Y sinc((x-shiftX(l))./pi)];
            end

            % add noise
            E = mvnrnd(zeros(1, V(j)), Omega, N(k));
            T = Y + E;

            Phi	= sbl_kernelFunction(x, x, kernelType, kernelWidth);

            % perform original MRVR or fast MRVR
            for l = 1:2
                tic
                switch l
                    case 1
                        % existing method
                        [used, ~, Mu, ~, OmegaHat, nIters(i,j,k,l)] ...
                            = mrvr(Phi,T,maxIts,tolerance);
                    case 2
                        % proposed method
                        [used, ~, Mu, ~, OmegaHat, nIters(i,j,k,l)] ...
                            = fmrvr(Phi,T,maxIts,tolerance);
                    otherwise
                        error('unsupported algorithm')
                end
                runTime(i,j,k,l) = toc;
                nRvs(i,j,k,l) = length(used); % # of relevance vectors

                Ypredict = Phi(:,used)*Mu;
                rmse(i,j,k,l) = sqrt(mean(mean((Y - Ypredict).^2)));

                % loss functions to measure the performance of covariance matrix estimation
                temp = (Omega\OmegaHat)'; % Omega is true, and OmegaHat is estimated
                loss1(i,j,k,l) = trace(temp) - log(det(temp)) - V(j); % entropy loss
                loss2(i,j,k,l) = trace((temp - eye(V(j)))^2); % quadratic loss
            end
        end
    end
end

%% median values
medRmse = squeeze(median(rmse));
medRunTime = squeeze(median(runTime));
medNRvs = squeeze(median(nRvs));
medLoss1 = squeeze(median(loss1));
medLoss2 = squeeze(median(loss2));
medNIters = squeeze(median(nIters));

%% mesh and table
close all
az = -37.5; % angle of view
el = 30; % angle of view
lineWidth = 2;

% row names and variable names of table
rowNames = mat2cell([repmat('V=',length(V),1) num2str(V')],...
    ones(1,length(V)),2+ceil(log10(max(V))));
variableNames = mat2cell([repmat('N',length(N),1) num2str(N')],...
    ones(1,length(N)),1+ceil(log10(max(N))));

% remove space in variableNames
for i = 1:length(N)
    temp = variableNames{i};
    variableNames{i} = temp(~isspace(temp));
end

% loop to generate mesh and table
for i = 1:6
    switch i
        case 1
            Z = nIters;
            medZ = medNIters;
            zScale = 'linear';
            titleName = 'The number of iterations';        
        case 2
            Z = runTime;
            medZ = medRunTime;
            zScale = 'log';
            titleName = 'Running time (seconds)';
        case 3
            Z = loss1;
            medZ = medLoss1;
            zScale = 'log';
            titleName = 'Entropy loss';
        case 4
            Z = loss2;
            medZ = medLoss2;
            zScale = 'log';
            titleName = 'Quadratic loss';
        case 5
            Z = rmse;
            medZ = medRmse;
            zScale = 'linear';
            titleName = 'RMSE';            
        case 6
            Z = nRvs;
            medZ = medNRvs;
            zScale = 'linear';
            titleName = 'The number of RVs';
    end
    figure('Position',[100 50 400 400]); % [left, bottom, width, height]
    h(1) = mesh(N,V,medZ(:,:,1),'linestyle','--'); hold on
    set(h(1),'facecolor','none')
    set(h(1),'LineWidth',lineWidth)
    h(2) = mesh(N,V,medZ(:,:,2)); hold off
    set(h(2),'facecolor','none') 
    set(h(2),'LineWidth',lineWidth)
    xlim([N(1) N(end)])
    xlabel('N')
    set(get(gca, 'xlabel'), 'FontAngle', 'italic')
    ylabel('V')
    set(get(gca, 'ylabel'), 'FontAngle', 'italic')
    set(gca, 'XTick', N, 'YTick', V)
    view(az,el)
    set(gca, 'ZScale', zScale)
    zlim([min(medZ(:)) max(medZ(:))])
    zlabel(titleName)
    set(h(1),'CData',zeros(length(V),length(N)));
    set(h(2),'CData',ones(length(V),length(N)));
    switch i
        case {3,4,5}
            hLeg = legend('Existing method','Proposed method');
            set(hLeg,'Location','NorthEast')
    end
    
    % difference in median values
    diff = medZ(:,:,1) - medZ(:,:,2);
    tableDiff{i} = array2table(diff,'RowNames',rowNames,...
        'VariableNames',variableNames);
    
    % p-value
    pval = nan(length(V), length(N));
    h = pval;
    for j = 1:length(V)
        for k = 1:length(N)
            [pval(j,k), h(j,k)] = ranksum(Z(:,j,k,1), Z(:,j,k,2));
        end
    end
    
    % table of 1) p-value and 2) rejection of null hypothesis
    tablePval{i} = array2table(pval,'RowNames',rowNames,...
        'VariableNames',variableNames);
    tableH{i} = array2table(h,'RowNames',rowNames,...
        'VariableNames',variableNames);
    
    % check normality using JB test
    nGrids = length(V).*length(N);
    for l = 1:2
        for j = 1:length(V)
            for k = 1:length(N)
                h(j,k) = jbtest(Z(:,j,k,l));
            end
        end
        tableJb{l,i} = array2table(h,'RowNames',rowNames,...
        'VariableNames',variableNames);
    
        nRejectNull(l,i) = sum(h(:));
    end
end

%% plot of true functions
figure, plot(x,Y,'LineWidth',lineWidth)
xlabel('x'), set(get(gca, 'xlabel'), 'FontAngle', 'italic', 'FontName', 'Times')
ylabel('y'), set(get(gca, 'ylabel'), 'FontAngle', 'italic', 'FontName', 'Times')
ylim([-.4 1.2])
legend('1st output','2nd output','3rd output','4th output','5th output')