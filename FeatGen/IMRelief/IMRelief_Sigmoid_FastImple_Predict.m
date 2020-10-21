function [pred, max_prob, max_ind] = ...
    IMRelief_Sigmoid_FastImple_Predict(tr_patterns, tr_targets, ts_patterns, Weight, distance, sigma)
%--------------------------------------------------------------------------
Uc = unique(tr_targets);
if min(Uc)==-1, tr_targets(tr_targets==-1) = 2; end

lUc = length(Uc); N = zeros(lUc,1); X = cell(lUc,1);
for c = 1:lUc,
    indC = tr_targets == c;
    N(c) = sum(indC); 
    X{c} = tr_patterns(:,indC); 
end 

[dim, N_ts_patterns] = size(ts_patterns);

Prediction = zeros(N_ts_patterns,1);
max_prob = zeros(N_ts_patterns, lUc);
max_ind = zeros(N_ts_patterns, lUc);

% Compute Probabilities
for i = 1:N_ts_patterns,
        
        
        for c = 1:lUc
    
            switch lower(distance)
                case {'euclidean'}
                    Temp = (X{c} - repmat(ts_patterns(:,i),1,N(c))).^2;
                case {'block'}
                    Temp = abs(X{c} - repmat(ts_patterns(:,i),1,N(c)));
            end

            dist    = (Weight(:).^2)'*Temp;
            temp_index = find(dist==0);
            
            %calculate probabilities
            prob = exp(-dist/sigma);
            
            %a = ( ( Weight(:).^2 )' ) * Temp; % compute induced margin
            %prob = 1 ./ ( 1 + exp(-a) ); % Probability of being class 1
            
            if ~isempty(temp_index);prob(temp_index(1)) = 0;end

            if sum(prob)~=0;
                prob_1 = prob/sum(prob);
            else
                [dum,I] = sort(dist);
                prob(I(2))=1;
                prob_1=prob;
            end
            
            [max_prob(i,c), max_ind(i,c)] = max(prob_1);
            
        end
        
end

[dum, pred] = max(max_prob,[],2);

if min(Uc)==-1
    pred(pred==2) = -1;
end