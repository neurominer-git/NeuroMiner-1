function [y] = ml_visualize_tSNE(X, options)
% Description:
%     Implements symmetric tSNE from van der Maatens & Hinton 2008. Draws 
%     on parameter initialization constants in paper, and speedup tricks 
%     from van der Maatens Matlab implementation.
%
% Options:
%     gradient descent tuning parameters
%     - initial_momentum: starting momentum rate for gradient descent 
%                         iterations. Default: 0.5
%     - final_momentum: final momentum rate for gradient descent iters.
%                        Default: 0.8
%     - momentum_switch_iter: which iteration to switch between initial
%                              and final momentum rates. Default: 250
%     - eta: step size for gradient descent. Default: 100;
%     - exaggerate_factor: early iterations artificially inflate the
%     magnitude of similarity between high dimensional points to increase
%     the. Default: 4;
%     - last_exaggeration_iter : Default: 50;
%     - descent_iters : Default: 1000;
%     - gamma : Default: .2; % for Delta-Bar-Delta Learning Rule (Jacobs 1988)
%     - kappa : Default: 0.2;
%     % cap the minimum reduction in learning rate
%     - min_rate : Default: .01;

if isempty(options)
    options = [];
end

[max_dims, tSNEOptions] = myProcessOptions(options, 'max_dims',50, ...
                                           'tSNEOptions', []);
                                       
% Pre-process model matrix using PCA to reduce dims to at most max_dims
X = center_and_reduce(X, max_dims);

% Find pairwise distances
D = squared_L2_norm(X);

% Find conditional probabilities
P = find_conditional_prob(D, options);

% Find lower dimensional representation of data
y = tSNE(P, tSNEOptions);

end

function [ydata] = tSNE(P, optOptions)
% Symmetric version of t-Distributed Stochastic Neighbour Embedding

% Relies on implementation tricks from paper and from Van der Maatens
% implementation. The goal of this re-implementation is to integrate
% t-SNE into the matLearn API for easy experimentation and comparison with 
% other ML approaches.

% See documentation of ml_visualize_tSNE for options documentation

nInstance = size(P, 1); % number of instances

[initial_momentum, final_momentum, momentum_switch_iter, eta, ...
 exaggerate_factor, last_exaggeration_iter, iters, gamma, kappa, ... 
 min_rate, verbose, final_dims] = myProcessOptions(optOptions, ...
                                         'initial_momentum', 0.5, ...
                                         'final_momentum', 0.8, ...
                                         'momentum_switch_iter', 250, ...
                                         'eta', 100, ...
                                         'exaggerate_factor', 4, ...
                                         'last_exaggeration_iter', 100, ...
                                         'descent_iters', 1000, ...
                                         'gamma', 0.2, ...
                                         'kappa', 0.2, ...
                                         'min_rate', 0.1, ...
                                         'verbose', 1, ...
                                         'final_dims', 2);

moment = initial_momentum;

% make conditional probability distribution symmetric for this algorithm
P = 0.5 * (P + P');

% guarantee sum-to-1 holds for probabilities after rounding errors
P = max(P ./ sum(P(:)), realmin); % van der Maatens trick                    

% "early exaggeration" to encourage optimization to focus on large p_{ij}
% by selecting large q_{ij}
P = P * exaggerate_factor;

% Sample initial solution in low-d space
ydata = 1e-4 * randn(nInstance, final_dims); % random sample from N(0,1e-4)

% preallocate data structures
[n_y, d_y] = size(ydata);
weights = ones(n_y, d_y);
y_step = zeros(n_y, d_y);

for iter=1:iters
    % Find gradient
    D_y = squared_L2_norm(ydata);
    t_density = 1./(1+D_y); % t-dist with 1 df, a.k.a. Cauchy distribution
    Q = t_density ./ sum(t_density(:)); % compute low-d affinities
    
    % fast gradient calculation from van der Maatens implementation of
    % t-SNE
    L = (P - Q) .* t_density;
    grad = 4 * (L - diag(sum(L, 1))) * ydata;

    % adapt gradient step for each point using Jones (1988)
    % delta-bar-delta rule
    weights = adapt_learning_rates(weights, grad, y_step, kappa, ...
                                   gamma, min_rate);
                               
    % Calculate adapted descent step down the gradient                           
    y_step = moment * y_step + eta * (grad .* weights);
    % Take step
    ydata = ydata + y_step;
    
    % recenter ydata
    ydata = ydata - repmat(mean(ydata),n_y,1); % van der Maatens trick
    
    if iter == last_exaggeration_iter
        P = P ./ exaggerate_factor;
    end
    
    if iter == momentum_switch_iter
        moment = final_momentum;
    end
end
end

function [X] = center_and_reduce(X, max_dims)
    [~, nVars] = size(X);
    X = standardizeCols(X);
    options_PCA = [];
    % Reduce high-dimensional space to at most max_dims dimensions
    options_PCA.maxComponents = min(max_dims, nVars);
    model_X = ml_unsupervised_dimRedPCA(X, options_PCA);
    X = model_X.reduceFunc(model_X, X);
    X = standardizeCols(X);
end

function [P] = find_conditional_prob(D, options)
[perp, maxIter, tol] = myProcessOptions(options, 'perplexity', 25, ...
                                        'tol', 1e-3, 'maxIter',50);
[nInstances, ~] = size(D);

% Do bracketed bisection search: finds precision tau_i acheiving target 
% entropy for induced probability distributions
targetEntropy = log(perp);
tau = ones(nInstances, 1);
P = zeros(nInstances);

for i = 1:nInstances
    % dists from high-dimensional x_i to x_j s.t. j != i
    D_i = D(i, [1:i - 1, i+1:end]);
    tau(i) = bisection_search(@compute_H_i, D_i, targetEntropy, maxIter, tol);
    P(i, [1:i-1, i+1:end]) = compute_P_i(D_i, tau(i)); % P_{i,i} := 0
end

end

function [rates] = adapt_learning_rates(rates, grad, y_step, kappa, ...
                                        gamma, min_rate)
    % Use delta-bar-delta update rule for adaptive learning rates: 
    % If continuing in the same direction at this step as last,
    % increment weight by kappa. If not, scale by 1 - gamma
    
    % Note: method calculated negative of gradient
    same = sign(grad) == sign(y_step);
    rates = (rates + kappa) .* same + (rates * (1 - gamma)) .* (~same);
    % enforce minimum rate constraint
    rates(rates < min_rate) = min_rate;
end

function [theta] = bisection_search(func, fixed_param, target, maxIter, tol)
    % Finds theta s.t. abs(f(theta)-target) < tol for a bivariate function
    % func(w,x) by holding w=fixed_param constant and searching for best x.
    % Requires func(w,x) to be monotone decreasing in x with w fixed
    
    iter = 0;
    theta = 1;
    % Initialize bracketing
    LB = -Inf;
    UB = +Inf;
    currentVal = func(fixed_param, theta);
    delta = target - currentVal;
    while (iter < maxIter)
        if delta < 0
            LB = theta;
            if isinf(UB)
                theta = 2*theta;
            else
                theta = (theta + UB)/2;
            end
        else
            UB = theta;
            if isinf(LB)
                theta = theta/2;
            else
                theta = (theta + LB)/2;
            end
        end
        currentVal = func(fixed_param, theta);
        delta = target - currentVal;
        if abs(delta) < tol
            break;
        end
        iter = iter + 1;
    end
end   

function [D] = squared_L2_norm(X)
    [n,d] = size(X);    
    D = X.^2*ones(d,n) + ones(n,d)*(X').^2 - 2*(X*X');
end

function [H] = compute_H_i(tau, D_i)
    % Compute the log perplexity (a.k.a. Shannon entropy) of the Gaussian 
    % distribution induced by precision paramater tau
    P_numer = exp(-D_i * tau / 2); % exp(-||x_i - x_j||^2 * tau / 2)
    P_denom = sum(P_numer); % sum_{k not i} exp(-||x_i - x_k||^2 * tau / 2)
    
    % form vector of exp(-||x_i - x_j||^2 * tau / 2) / P_denom for each j
    P_all = D_i ./ P_denom;
    
    % calculate Shannon entropy using a simplified form of the equation by
    % subsituting natural log for binary log
    H = log(P_denom) + tau * sum(P_all) / P_denom;
end

function [P_i] = compute_P_i(D_i, tau)
% Compute similarity of datapoints x_j to x_i as Pj|i given precision 
% tau (inverse variance) for each j != i 
    P_i = exp(-D_i * tau / 2) / sum(exp(-D_i * tau / 2));
end
