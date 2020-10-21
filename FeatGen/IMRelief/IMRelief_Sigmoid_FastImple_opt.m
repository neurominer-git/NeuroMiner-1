function [Weight, tElapsed] = IMRelief_Sigmoid_FastImple(patterns, targets, distance, sigma, lambda, tmax, plotfigure)
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
%Modified by Nikolaos Koutsouleris, 06/2011

fprintf('\nParameters: distance metric = %s, sigma = %g, lambda = %g', distance, sigma, lambda);
%--------------------------------------------------------------------------
Uc = unique(targets);
if min(Uc)==-1
    %targets(targets==1) = 2;
    %targets(targets==-1) = 1;
    targets = targets/2+1.5; % transform targets into {1,2}
end
Targets = targets-1;
[dim,N_patterns] = size(patterns);      % Data dimenionality

Original_dim = dim;
Original_index = 1:dim;
index_0 = find(Targets==0);
History = zeros(dim,tmax+1);
Weight = ones(dim,1);
History(:,1) = Weight;
T = Targets(:);
% -------------------------------------------------------------------------

Difference = 1;
t=0;
theta = zeros(tmax,1); tic

% Define X
lUc = length(Uc); N = zeros(lUc,1);
X = cell(lUc,1); Y = cell(lUc, N_patterns);
for c = 1:lUc,
    indC = targets == c;
    N(c) = sum(indC); 
    X{c} = patterns(:,indC); 
     for n = 1:N_patterns
         Y{c,n} = repmat(patterns(:,n),1,N(c));
     end
end 

while  Difference>0.01 && t<=tmax;
    
    if t>0, fprintf('\nIter #%g: theta = %1.2f',t, Difference);end
    t=t+1;
    NM = zeros(dim, N_patterns);
    NH = zeros(dim, N_patterns);
    if t>1, fprintf(' => Compute NM and NH ...'); end
   
    for i = 1:N_patterns,

        o=1; 
        
        for c = 1:lUc
           
            switch lower(distance)
                case {'euclidean'}
                    Temp = (X{c} - Y{c,i}).^2;
                    %Temp = (X{c} - repmat(patterns(:,i),1,N(c))).^2;
                case {'block'}
                    Temp = abs(X{c} - Y{c,i});
                    %Temp = abs(X{c} - repmat(patterns(:,i),1,N(c)));
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
                    mTemp = Temp*prob_1(:); 
                    mdist = min(dist); 
                end
                o=o+1;
            end
        end
        
        NM(:,i) = NM(:,i) + mTemp;
    end
    
    Z = NM-NH;
    Z(:,index_0) = -Z(:,index_0);
    Cost = zeros(1,2); CostDiff = 1000; Cost(2) = 10000; 
    if t > 1, fprintf(' => Update feature weights ...'), end
    
    while CostDiff>0.01*Cost(2) 
        
        a = ((Weight(:).^2)')*Z; % Margin
        Result = 1./(1+exp(-a)); % Probability of being class 1
        difference = Result(:)-Targets(:);
        descent = (Z*difference(:)).*Weight+lambda*Weight;

        %line search to find alpha    
        Cost(1) = Cost(2);
        [alpha, Cost(2)] = ...
            fminbnd('IMRelief_cost', 0, 1, optimset('TolX',0.02), Z, Weight, descent, T, lambda); 

        Weight = Weight-alpha*descent;
        CostDiff = abs(Cost(2)-Cost(1));
       
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
             for n = 1:N_patterns
                 Y{c,n} = Y{c,n}(indkeep,:);
             end
        end
        patterns = patterns(indkeep,:);
        dim = sum(indkeep);
        Weight = Weight(indkeep);
        History = History(indkeep,:);
        Original_index = Original_index(indkeep);
    end
    
    if t > 1, fprintf(' Done.'); end
end
tElapsed = toc;
fprintf('\nFinal theta = %1.2f. Done in %1.2f sec.',Difference, tElapsed) 
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

return

