function [Objective, Weight, History] = iDetect(patterns, Para)
% iDetect: feature selection for clustering and manifold learning with 
%           norm-2 and norm-1 regularization
%INPUT:
%        patterns:  D-by-N input data
%                   D is # of features and N is # of samples 
%            Para:  parameter structure
%         Para.it:  maximum iteration     
%   Para.distance:  distance metric. 'euclidean' or 'block' (default)
%      Para.sigma:  kernel width
%     Para.lambda:  regularization parameter
%                   (can be tuned by gap statistics, see reference)
%                   
% 
%OUTPUT:
%          Weight: weight of features 
%       Objective: objective function value
%         History: feature weight in every iteration
%
%REFERENCE:
% J. Yao, Q. Mao, V. Mai, S. Goodison, Y. Sun, 
% Feature Selection for Unsupervised Learning through Local Learning
% Pattern Recognition Letters, 2015
%
% VERSION: 12/22/2014
%--------------------------------------------------------------------------
%% ==========================================================================
if(nargin==1)
    T = 20;                            
    distance = 'block'; 
    sigma = 1;       
    lambda = 2;
else
    T = Para.maxiter;                            
    distance = Para.distance; 
    sigma = Para.sigma;       
    lambda = Para.lambda;     
end

%% normalize the data so that each feature is comparable
[dim,N_patterns] = size(patterns);
% Normalize to [0, 1]
[MIN,~] = min(patterns,[],2);
[MAX,~] = max(patterns,[],2);
for n=1:dim
    if(MIN(n)==MAX(n))
        patterns(n,:) = 0;
    else
        patterns(n,:) = (patterns(n,:)-MIN(n))/(MAX(n)-MIN(n));
    end
end

%% initilization
Weight = 1/sqrt(dim)*ones(dim,1); %initial guess
History = zeros(dim,T+1);
History(:,1) = Weight;
Objective = [];
Difference = 1;
t = 0;
Theta = [];

Original_dim = dim;
Original_index = 1:dim;
OriginalPatterns = patterns;
%% main body of the code
while Difference>0.01 && t<=T; 
    t=t+1;
    NN = zeros(dim,N_patterns);
    FN = zeros(dim,N_patterns);
    
        for i = 1:N_patterns,
            switch lower(distance)
                case {'euclidean'}
                    Temp            = (patterns - patterns(:,i)*ones(1,N_patterns)).^2;
                case {'block'}
                    Temp            = abs(patterns - patterns(:,i)*ones(1,N_patterns));
            end

             dist    = (Weight(:)')*Temp;
            [~,temp_index] = sort(dist);

            % probability of being the nearest neighbor
            prob = exp(-(dist-dist(temp_index(2)))/sigma);prob(temp_index(1)) = 0;
            prob_1 = prob/sum(prob);
            NN(:,i) = Temp*prob_1(:);
            FN(:,i) = sum(Temp,2)/(N_patterns-1);
        end
        
    C = sum(FN-NN,2);C_original = C;
    C(C<0) = 0;
    Weight = C/norm(C);
    
    if sum(Weight)>lambda
        min_delta = 0;
        max_delta = max(C_original);
        min_delta_dif = 10^(-3); % resolution, stop if the difference is less than this number.
        min_conv_dif = 10^(-3);  % convergence stopping, if abs(weightsum-old_weightsum) less than this number
        dif_conv = 1;%sum(Weight)-lambda; %to verify convergence stopping
        dif_delta = 1;%max_delta-min_delta;
        while dif_delta > min_delta_dif && dif_conv > min_conv_dif
            mid_delta = min_delta + (max_delta - min_delta)/2;          
            C = C_original - mid_delta;
            C(C<0) = 0;
            Weight = C/norm(C);
            if sum(Weight) > lambda
                min_delta = mid_delta;
            end
            if sum(Weight) <= lambda
                max_delta = mid_delta;
            end
            dif_conv = abs(sum(Weight)-lambda);
            dif_delta = max_delta - min_delta;
        end
    end
      
    Objective = [Objective, (Weight(:)')*C_original];
    
    temp = zeros(Original_dim,1);
    temp(Original_index) = Weight;
    Weight = temp;
    Original_index = 1:Original_dim;

    Difference = norm(Weight-History(:,t));
    Theta(t) = Difference;
    History(:,t+1) = Weight;
       
    %a huristic to improve the speed 
    if t>4
        SUM = sum(History(:,t-4:end),2);
        Zero_index = find(SUM==0);
        patterns = OriginalPatterns;
        patterns(Zero_index,:) = []; %if the sum of a weight in the past three iteration is zero, the feature is removed.
        Original_index(Zero_index) = [];
        dim = size(patterns,1);
        Weight(Zero_index) = [];
    end
   
end

% Plot objective history
% figure;plot(Objective,'-o')
% title('Objective function');
% xlabel('Number of Iteration');
% ylabel('Objective function')
% grid on
% hold on

Weight = History(:,t+1);
History = History(:, 1:t+1);
Objective = Objective(end);

% monitoring the feature weights
% figure;plot(Theta,'-o')
% title('Theta');
% xlabel('Number of Iteration');
% ylabel('Difference')
% grid on
% hold on
% figure;semilogx(Weight,'-ro')
