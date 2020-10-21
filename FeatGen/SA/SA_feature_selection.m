function SA_feature_selection(Y, label, Ynew, labelnew)

% #####################################################################
% This code does:
% 1) Feature selection using SA
% 2) This is a simple code working on 1 single shot only. That is, this
% code work with a set of trainin, cross-validation and test set. 
% 3) If you want to make N-fold cross validation, you can certainly do so
% by covering the code by a for-loop.
% #####################################################################

 n = size(Y,2); % the number of data instances (m) and the features (n)

% #####################################################################
% ==== Make the order of voxels to be introduced to the SA =====
% #####################################################################
% The normal order
featureID_sorted = 1:n;
% featureID_order = featureID_sorted(randperm(n)); % Randomly shuffle the voxel order
featureID_order = featureID_sorted ; % the MI-descend order
%% 
% =================================================================
%                       SA framework
% =================================================================
% ===== SA operating parameters ==== @#$% user-defined
% @#$% user-defined objective function
objFunction = @ObjectiveFunction1; 

% @#$% user-defined solution mapping function
nextSolution = @NextSolution3; 

% @#$% SA parameters
T               = 1;
T_stop          = 0.005; % the stopping temperature 
alpha           = 0.9;
itt_max         = 2000;
Rep_T_max       = 100; % max number of Rep in one temperature T
Rep_accept_max  = 30; % if solution accepted this many, then update T
kc              = 10; % k-constant. The smaller kc --> less solution accepted @#$% user-defined

% ====================================================
% initialize or evaluate the solution for user-defined function  % @#$% user-defined
x_curr          = zeros(1,n); x_curr(1:5) = 1; % starting with choosing first 5 voxels
F_curr          = -1e+9; % F_curr = objFunction(x_curr);
x_best          = x_curr;
F_best          = F_curr;
storage_new     = {};
% ====================================================

% ====================================================
% initial value for SA
Rep_T           = 0;
Rep_accept      = 0;
itt             = 1;
converges       = 0;
is_best_updated = 0; % toggle flag set to 1 when the best  solution is updated

% Start SA optimization
while converges == 0 && itt <= itt_max
    
    % pick a new solution x_new
    x_new = nextSolution(x_curr, itt, n, featureID_order); % @#$% user-defined
    
    % @#$% user-defined objective function
    % F_new = objFunction(x_new): Calculate F_new from x_new: 
    F_new = objFunction(n, c, labelnew, Ynew(:,x_new), label, Y(:,x_new));
    
    % Anything you want to keep is here % @#$% user-defined
    storage_new.itt = itt;
    storage_new.sol = x_new;
    storage_new.value = F_new;
    storage_new.alpha = alpha;
    storage_new.T = T;
        
    % The code here is exactly the same as the one above except not
    % showing any plot which is much faster.
    % Check: update solution
    if F_new >= F_curr
        % back up before update
        F_prev = F_curr; % backup the previous
        x_prev = x_curr; % backup the previous
        % update
        F_curr = F_new; % accept the solution
        x_curr = x_new; % accept the solution
        Rep_accept = Rep_accept + 1;
        % check the best solution
        if F_new >= F_best
            F_best = F_new;
            x_best = x_new;
            storage_best = storage_new;
            is_best_updated = 1;
        end
    elseif exp((F_new-F_curr)/(kc*T)) > rand(1)
        % back up before update
        F_prev = F_curr; % backup the previous
        x_prev = x_curr; % backup the previous
        % update
        F_curr = F_new; % accept the solution
        x_curr = x_new; % accept the solution
        Rep_accept = Rep_accept + 1;
    end

    % update the counters
    Rep_T = Rep_T + 1;
    itt = itt + 1;
    
    % Check: update temperature
    if Rep_T >= Rep_T_max || Rep_accept >= Rep_accept_max
        T = alpha*T;
        Rep_T = 0;
        Rep_accept = 0;
        F_curr = F_best;
        x_curr = x_best;
    end
    
    % stop criteria by temperature
    if T < T_stop
        converges = 1;
    end
    
    
    % %     % Check for convergence???
    % %     % -- Well, in some objective function, it's very hard to define the
    % %     % stopping criteria because the solution is very jumy even at the end,
    % %     % therefore, I don't check convergence here. Anyway, you can do so
    % %     % here.
    % %     % =========================================================
    % %     if F_curr - F_prev > something
    % %         converges = 1;
    % %     end
    
end

% =========== END of SA =====================
toc;