function [Pred, Pro] = ...
    IMRelief_Sigmoid_FastImple_Predict2(train_patterns, train_targets, test_patterns, Weight, distance, sigma)

index_1 = find(train_targets==-1);N(1) = length(index_1);
index_2 = find(train_targets==+1);N(2) = length(index_2);
patterns_1 = train_patterns(:,index_1);
patterns_2 = train_patterns(:,index_2);
Weight = Weight(:);
[dim, N_test_patterns] = size(test_patterns);
Pro = zeros(N_test_patterns,1);

for n = 1:N_test_patterns
    
    test = test_patterns(:,n);
    
    switch distance
        case 'euclidean'
             temp_1 = ( patterns_1 - test * ones( 1, N(1) ) ).^2;
             temp_2 = ( patterns_2 - test * ones( 1, N(2) ) ).^2;
        case 'block'
             temp_1 = abs( patterns_1 - test * ones( 1, N(1) ) );
             temp_2 = abs( patterns_2 - test * ones( 1 ,N(2) ) );
    end
   
    dist_1 = (Weight)'*temp_1;
    prob_1 = exp(-dist_1/sigma);
    prob_1 = prob_1/sum(prob_1);
   
    
    dist_2 = (Weight)'*temp_2;
    prob_2 = exp(-dist_2/sigma);
    prob_2 = prob_2/sum(prob_2);
   
    rho = sum(dist_1.*prob_1)-sum(dist_2.*prob_2);
    Pro(n) = 1/(1+exp(-rho));
end

Pred = zeros(size(Pro));
Pred(Pro>=0.5) = 1; Pred(Pro<0.5) = -1;