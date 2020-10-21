function [model] = ml_unsupervised_LGSSM(X, options)
% Linear-Gaussian State Space Models
%
%   Implements Kalman filtering and Kalman smoothing for directed
%   probabilistic models that have Gaussian latent variables with a Markov
%   transition structure and Gaussian observation structure.

[A, C, mu_init, V_init, G, S] = myProcessOptions(options, 'A', [], ...
                                                 'C', [], ...
                                                 'mu_init', [], ...
                                                 'V_init', [], 'G', [], ...
                                                 'S', []);

model.KalmanFilter = @(x)KalmanFilter(A, C, G, S, x, mu_init, V_init);
model.KalmanSmoothing = @(x, y, z) KalmanSmoothing(A, G, C, x, y, z);
end

function [mu, V, K] = KalmanFilter(A, C, Q, S, X, mu_0, V_0) 
[latentDim] = size(A, 1);
[seqLen] = size(X, 2);
mu = zeros(latentDim, seqLen);
V = zeros(latentDim, latentDim, seqLen);
for t = 1:seqLen
    if (t == 1) % start with initialization mean
        mu_pred = mu_0;
        V_pred = V_0;
    else % form next mean as linear combination of previous means
        mu_pred = A*mu(:, t-1); % mu_pred = x_t^T, for E[x_t | {y}_1^tau]
        V_pred = A*V(:, :, t-1)*A' + Q;
    end
    K = V_pred*C' / (C*V_pred*C' + S); % Kalman gain
    mu(:, t) = mu_pred + K*(X(:, t) - C*mu_pred); % correct prediction
    V(:, :, t) = V_pred - K*C*V_pred;
end
end

function [mu_hat, V_hat, P_t, P_t_t1] = KalmanSmoothing(A, G, C, mu, V,K_T)
[latentDim] = size(A, 1);
[seqLen] = size(mu, 2);
mu_hat = zeros(latentDim, seqLen);
V_hat = zeros(latentDim, latentDim, seqLen);
V_hat_T = zeros(latentDim, latentDim, seqLen);
J = zeros(latentDim, latentDim, seqLen);
P_t = zeros(latentDim, latentDim, seqLen);
P_t_t1 = zeros(latentDim, latentDim, seqLen); 


mu_hat(:, seqLen) = mu(:, seqLen);
V_hat(:,:,seqLen) = V(:,:,seqLen);
% Implements Hinton and Ghahramani (G & H) 1996, Parameter Estimation of 
% Linear Dynamical Systems
V_hat_T(:, :, seqLen) = (eye(latentDim) - K_T * C) * A * ...
                                                     V(: ,:, seqLen - 1);

for t = seqLen - 1:-1:1
    % Useful intermediate from Forward pass:
    M_t = A*V(:, :, t)*A' + G;
    J(:, :, t) = V(:, :, t)*A' / M_t; 
    mu_hat(:, t) = mu(:, t) + J(:,:,t)*(mu_hat(:, t+1) - A*mu(:, t));
    V_hat(:, :, t) = V(:, :, t) + J(:,:,t)*(V_hat(:,:,t+1)- M_t)*J(:,:,t)';
end
% Calculate a useful intermediate for EM iterations, based on G & H 1996
P_t(:,:,seqLen) = V_hat(:,:,seqLen) + mu_hat(:,seqLen)*mu_hat(:,seqLen)';
P_t_t1(:,:,seqLen) = V_hat_T(:,:,seqLen) + mu_hat(:, seqLen)* ...
                                           mu_hat(:, seqLen - 1)';
for t = seqLen - 1: -1 : 2 
    V_hat_T(:,:,t) = V(:,:,t)*J(:,:,t-1)' + J(:,:,t) * ... 
                     (V_hat_T(:,:,t+1) - A*V(:,:,t))*J(:,:,t-1)';
    P_t(:,:,t) = V_hat(:,:,seqLen) + mu_hat(:,t)*mu_hat(:,t)';
    P_t_t1(:,:,t) = V_hat_T(:,:,t) + mu_hat(:,t)*mu_hat(:,t - 1)';
end
P_t(:,:,1) = V_hat(:,:,1) + mu_hat(:,1)*mu_hat(:,1)';
end

