function [u, v, success] = spls(X, Y, cu, cv, e, itr_lim)
%
%   Sparse PLS algorithm, please check Monteiro et al. 2016 for details:
%   doi:10.1016/j.jneumeth.2016.06.011
%
%   Inputs: X, Y    - data matrices in the form: samples x features. These
%                     should have each feature with mean = 0 and std = 1;
%
%           cu, cv  - sparcity regularization hyperparameters, must be
%                     between 1 and sqrt(number_features). The lower it is,
%                     the spaser the solution. If it is outside this range,
%                     no sparsity will be applied in the corresponding view.
%
%           e       - convergence threshold (see the code for info on how it
%                     works). Default: 1E-5
%
%           itr_lim - maximum number of iterations (it give a warning
%                     if it does not converge). Default: 1000
%
%
%   Outputs: u, v    - weight vectors for X and Y, respectively
%
%            success - will return "false" if something went wrong during
%                      the weight vector computation
%
%   Version: 2016-08-20
%__________________________________________________________________________

% Written by Joao Matos Monteiro
% Email: joao.monteiro@ucl.ac.uk


%--- Initial checks
%--------------------------------------------------------------------------

% Check if lu anv lv obey the limits
global VERBOSE
if cu < 1 || cu > sqrt(size(X,2))
    warning('lu is out of interval: 1 <= lu <= sqrt(size(X,2). Not using spasity on u.')
    no_sparse_X = true;
    failed_sparsity_u = false;
else
    no_sparse_X = false;
end
if cv < 1 || cv > sqrt(size(Y,2))
    warning('lv is out of interval: 1 <= lv <= sqrt(size(Y,2). Not using spasity on v.')
    no_sparse_Y = true;
    failed_sparsity_v = false;
else
    no_sparse_Y = false;
end

% Convergence threshold
if ~exist('e', 'var')
    e = 1E-5;
end

% Iteration limit for calculating a vector pair
if ~exist('itr_lim', 'var')
    itr_lim = 1000;
end



%--- SPLS
%--------------------------------------------------------------------------

%--- Compute the covariance matrix
C = X'*Y;


%--- Initialise weight vectors
u_temp = nan(size(X, 2), 2);
v_temp = nan(size(Y, 2), 2);


[U,~,V] = svd(C,0);
u_temp(:,1) = U(:,1);
u_temp(:,1) = u_temp(:,1)./norm(u_temp(:,1)); % normalise
v_temp(:,1) = V(:,1);
v_temp(:,1) = v_temp(:,1)./norm(v_temp(:,1)); % normalise

clear U V

%--- Main Loop
diff = 10*e; %start the diff with a high value
i = 0;
success = true;

while diff > e && success
    
    %--- Compute u
    if no_sparse_X
        u_temp(:,2) = C*v_temp(:,1);
        u_temp(:,2) = u_temp(:,2)./norm(u_temp(:,2), 2);
    else
        [u_temp(:,2), tmp_success] = update(C*v_temp(:,1), cu);
        failed_sparsity_u = ~tmp_success;
        if failed_sparsity_u % If it was not successful, return non sparse version
            u_temp(:,2) = C*v_temp(:,1);
            u_temp(:,2) = u_temp(:,2)./norm(u_temp(:,2), 2);
        end
    end
    dim_u = sum(u_temp(:,2)~=0);
    if ~dim_u
        error(['No weights were included in the model, this should never '...
            'happen. Try increasing lu.']);
    end
    
    
    %--- Compute v
    if no_sparse_Y
        v_temp(:,2) = C'*u_temp(:,2);
        v_temp(:,2) = v_temp(:,2)./norm(v_temp(:,2), 2);
    else
        [v_temp(:,2), tmp_success] = update(C'*u_temp(:,2), cv);
        failed_sparsity_v = ~tmp_success;
        if failed_sparsity_v % If it was not successful, return non sparse version
            v_temp(:,2) = C'*u_temp(:,2);
            v_temp(:,2) = v_temp(:,2)./norm(v_temp(:,2), 2);
        end
    end
    dim_v = sum(v_temp(:,2)~=0);
    if ~dim_v
        error(['No weights were included in the model, this should never '...
            'happen. Try increasing lv.']);
    end
    
    
    %--- Check convergence
    diff_u = norm(u_temp(:,2) - u_temp(:,1));
    diff_v = norm(v_temp(:,2) - v_temp(:,1));
    if diff_u >= diff_v, diff = diff_u; else diff = diff_v; end
    % update u and v for the next iteration
    u_temp(:,1) = u_temp(:,2);
    v_temp(:,1) = v_temp(:,2);
    
    if i >= itr_lim
        warning('Maximum number of iterations reached.');
        success = false;
    end
    
    i = i+1;
end

if failed_sparsity_u
    warning(['There was a problem with the delta estimation in u.' ...
        ' The solution was forced to be non-sparse. Take results with a grain of salt.']);
    success = false;
end

if failed_sparsity_v
    warning(['There was a problem with the delta estimation in v.' ...
        ' The solution was forced to be non-sparse. Take results with a grain of salt.']);
    success = false;
end

if VERBOSE; fprintf('\rSPLS: itr: %d    diff: %.2e    dim_u: %d    dim_v: %d', i, diff, dim_u, dim_v); end

%--- Add converged weight vectors to output
u = u_temp(:, end);
v = v_temp(:, end);


end


%--- Private functions
%--------------------------------------------------------------------------
function [up, success] = update(w, c)

success = true;

%--- update values
delta = 0;
up = soft_thresh(w, delta);
up = up./norm(up,2);

%--- check if it obeys the condition. If not, find delta that does.
if norm(up, 1) > c
    
    delta1 = delta;
    delta2  = delta1+1.1; % delta2 must be > 1
    
    % get first estimate of delta2
    flag = false;
    i = 0;
    max_delta = 0;
    while ~flag
        up = soft_thresh(w, delta2);
        up = up./norm(up,2);
        
        if sum(abs(up)) == 0 || isnan(sum(abs(up))) % if everthing is zero, the up/|up| will be 0/0 = nan
            delta2 = delta2/1.618; % They have to be diferent, otherwise it might not converge
        elseif norm(up, 1) > c
            delta1 = delta2;
            delta2 = delta2*2; % They have to be diferent, otherwise it might not converge
        elseif norm(up, 1) <= c
            flag = true;
        end
        
        if delta2>max_delta, max_delta = delta2;end
        
        if delta2 == 0
            warning('Delta has to be zero.');
            success = false;
            break
        end
        i = i+1;
        if i>1E4
            warning('First delta estimation update did not converge.');
            delta1 = 0;
            delta2 = max_delta;
            break
        end
    end
    

    up = bisec(w, c, delta1, delta2);
    if isempty(up) || sum(isnan(up))>0
        warning('Delta estimation unsuccessful.')
        success = false;
    end
    
    
end



end

function out = soft_thresh(a,delta)
% Performs soft threshold (it does not normalize the output)
diff = abs(a)-delta;
diff(diff<0) = 0;
out = sign(a).*diff;

end


function out = bisec(K, c, x1,x2)
converge = false;
success = true;
tolerance = 1E-6;
while ~converge && success
    x = (x2 + x1) / 2;
    out = soft_thresh(K, x);
    out = out./norm(out,2);
    if sum(abs(out)) == 0
        x2 = x;
    elseif norm(out, 1) > c
        x1 = x;
    elseif norm(out, 1) < c
        x2 = x;
    end
    
    diff = abs(norm(out, 1) - c);
    if diff <= tolerance
        converge = true;
    elseif isnan(sum(diff))
        success = false;
        out = nan(size(K));
    end
end
end