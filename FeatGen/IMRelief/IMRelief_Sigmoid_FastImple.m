function [Weight, tElapsed] = IMRelief_Sigmoid_FastImple(patterns, targets, distance, sigma, lambda, tmax, plotfigure, verbose)
%IMRelief_Sigmoid_FastImple: fast implemetation of IMRelief_Sigmoid
%--------------------------------------------------------------------------
%INPUT:
%     patterns:  training data: [x1,x2,...xn] Each column is an observation
%      targets:  class label = {1,2}
%     distance:  distance metric (default: block distance)
%        sigma:  kernel width
%       lambda:  regulariztion parameter
%   plotfigure:  1: plot of the learning process; 0: do not plot
%OUTPUT:
%       Weight:  weight of features
%--------------------------------------------------------------------------
%by Yijun Sun @University of Florida
%update history: Feb. 10/March 20, 2007/OCT 10, 2010
%==========================================================================
%Modified by Nikolaos Koutsouleris, 10/2011

if ~exist('verbose','var'), verbose = false; end
if verbose, fprintf('\nParameters: distance metric = %s, sigma = %g, lambda = %g', distance, sigma, lambda); end

%--------------------------------------------------------------------------
Uc = unique(targets);
if min(Uc)==-1
    targets(targets==1) = 1;
    targets(targets==-1) = 2;
end
Targets = targets-1; % 0 for class 1, 1 for class 2
Ucz = unique(Targets);
[dim,N_patterns] = size(patterns);      % Data dimenionality
Original_dim = dim;
Original_index = 1:dim;
%index_0 = find(Targets==0);
History = zeros(dim,tmax+1);
Weight = ones(dim,1);
History(:,1) = Weight;
T = Targets(:);
% -------------------------------------------------------------------------

Difference = 1;
t=0;
theta = zeros(tmax,1); 

% Define X
lUc = length(Uc); N = zeros(lUc,1); X = cell(lUc,1);
for c = 1:lUc,
    indC = targets == c;
    N(c) = sum(indC); 
    X{c} = patterns(:,indC); 
end 

while  Difference > 0.01 && t <= tmax;
    
    if t>0 && verbose, fprintf('\nIter #%g: theta = %1.2f',t, Difference);end
    t=t+1;
    NM = zeros(dim, N_patterns);
    NH = zeros(dim, N_patterns);
    if t>1 && verbose, fprintf(' => Compute NM and NH ...'); end
   
    for i = 1:N_patterns,

        o=1; 
        %fprintf('\n %g of %g',i,N_patterns);
        
        for c = 1:lUc
           
            switch lower(distance)
                case {'euclidean'}
                    %Temp = (X{c} - repmat(patterns(:,i),1,N(c))).^2;
                    Temp = bsxfun(@minus, X{c}, patterns(:,i)).^2;
                case {'block'}
                    %Temp = abs(X{c} - repmat(patterns(:,i),1,N(c)));
                    Temp = abs(bsxfun(@minus, X{c}, patterns(:,i)));
            end

            if t==1
                dist = sum(Temp,1)/sqrt(dim);
            else
                dist = (Weight(:).^2)'*Temp;
            end
            temp_index = find(dist==0);
            
            %calculate probabilities
            prob = exp(-dist/sigma);
            if ~isempty(temp_index);prob(temp_index(1)) = 0;end
            if sum(prob)~=0;
                prob_1 = prob/sum(prob);
            else
                [dum,I] = sort(dist);
                prob(I(2))=1;
                prob_1=prob;
            end
            
            % Nearest Hit
            if targets(i) == c;
                NH(:,i) = Temp * prob_1(:);
            % Nearest Miss
            else % Find min_(c=1...c) dist(x_n, NM^(c)(x_n(|w) across c groups
                if o==1 || min(dist) < mdist
                    mTemp = Temp * prob_1(:); 
                    mdist = min(dist); 
                end
                o=o+1;
            end
        end
        
        NM(:,i) = NM(:,i) + mTemp;
    end
    
    % Compute margin from nearest-miss and nearest-hit
    Z = NM-NH;
    
    % Check whether One-vs-All multi-group decomposition needs to be
    % performed:
    if numel(Ucz) == 2, 
        Ux = 1;
    else
        Ux = numel(Ucz);
    end
    
    if t > 1 && verbose, fprintf(' => Update feature weights ...'), end
    % Perform One-vs-ALL optimization if this is necessary (multi-group mode)

    for uz = 1:Ux
        
        index_0 = Targets == Ucz(uz);
        xTargets = zeros(size(Targets)); xTargets(index_0) = 0; xTargets(~index_0) = 1;
        uZ = Z; uZ(:,index_0) = -Z(:,index_0);
        
        if verbose && t > 1
            Zscore = sum(uZ(:));
            fprintf(' Margin = %g. ', abs(Zscore))
        end
        
        Cost = zeros(1,2); CostDiff = 1000; Cost(2) = 10000; 
        
        while CostDiff > 0.01*Cost(2) 

            a = ( ( Weight(:).^2 )' ) * uZ; % compute induced margin
            Result = 1 ./ ( 1 + exp(-a) ); % Probability of being class 1
            difference = Result(:) - xTargets(:);
            descent = (uZ * difference(:) ) .* Weight + lambda * Weight;

            %line search to find alpha    
            Cost(1) = Cost(2);
            [alpha, Cost(2)] = ...
                fminbnd('IMRelief_cost', 0, 1, optimset('TolX',0.02), uZ, Weight, descent, xTargets(:), lambda); 

            Weight = Weight - alpha * descent;
            CostDiff = abs(Cost(2)-Cost(1));

        end
    end
    
    Weight = abs(Weight);
    Difference = norm(abs(Weight/max(Weight)-History(:,t)/max(History(:,t))));
    theta(t) = Difference;
    History(:,t+1) = Weight;  
    
    % Eliminate features with negligible weight
    if t==1; index_zeros = Weight <= 10^(-5); end
    if t>=2; index_zeros = Weight <= 10^(-5); end
    
    if any(index_zeros)
        indkeep = ~index_zeros; 
        for c=1:lUc
            X{c} = X{c}(indkeep,:);
        end
        patterns = patterns(indkeep,:);
        dim = sum(indkeep);
        Weight = Weight(indkeep);
        History = History(indkeep,:);
        Original_index = Original_index(indkeep);
    end
    
    if t > 1 && verbose, fprintf(' Done.'); end
end

if isnan(Difference)
    warning('Theta is NaN. Check your data!')
    fprintf('Returning weight vector of the previous non-NaN iteration.')
    Weight = History(:,t);
end
if verbose,     
    fprintf('\nFinal theta = %1.2f.',Difference); 
    fprintf('\n# of remaining features: %g (%1.2f%%)\n', dim, dim*100/Original_dim);
end

temp = zeros(1,Original_dim);
temp(Original_index) = Weight.^2;
Weight = temp;

%Monitoring the feature weights
if plotfigure ==1
    figure;
    semilogy(theta,'-o','LineWidth',1,'MarkerFaceColor','w','MarkerSize',10)
    title('Theta');
    xlabel('Number of Iterations');
    ylabel('Difference')
    grid on
    boldify1
    drawnow
end


