% #####################################################################
% This code does:
% 1) Feature selection using SA
% 2) This is a simple code working on 1 single shot only. That is, this
% code work with a set of trainin, cross-validation and test set. 
% 3) If you want to make N-fold cross validation, you can certainly do so
% by covering the code by a for-loop.
% #####################################################################

%%
clear
close all
clc

%%
% #####################################################################
% Add necessary toolboxes
% #####################################################################
tic;
% addpath to the libsvm toolbox
addpath('../../../libsvm-3.12/matlab');
addpath('../../../use_libsvm');
% addpath to the miscellaneous
addpath('../../../toolbox_misc');
%%
% #####################################################################
% ===== Load the data and make some necessary info ======
% #####################################################################
dirData = '/share/Bot/Research/mvpa_data/matlab_format/gamma2p5_wholebrain';
dataFilename = 'Haxby8_betaMap_JoeFormat_vITmask';
load(fullfile(dirData,dataFilename));

clear betaMap;
clear classMatrix;

[m n] = size(data.beta); % the number of data instances (m) and the features (n)

% make run number for each instance
run = []; for i = 1:10, run = [run;i*ones(8,1)]; end

%% 
% #####################################################################
% Preprocess the data
% #####################################################################

% Normalize the feature
[X_norm, mu_X, sigma_X] = featureNormalize(data.beta);

% Do not need to sort the class labels
y = data.classLabel;
%% 
% #####################################################################
% Separate the data into train, validation and test set
% #####################################################################

% Get the test, validation (cv) and the training set
f = 1;
test_run_number = mod([f],10)+1;
cv_run_number = mod([f+1:f+3],10)+1;
train_run_number = setdiff([1:10],union(test_run_number,cv_run_number));
% The test set
idx_test = zeros(m,1);
for i = 1:length(test_run_number)
    idx_test = idx_test | (run == test_run_number(i));
end
% The cv set
idx_cv = zeros(m,1);
for i = 1:length(cv_run_number)
    idx_cv = idx_cv | (run == cv_run_number(i));
end
% The train set
idx_train = zeros(m,1);
for i = 1:length(train_run_number)
    idx_train = idx_train | (run == train_run_number(i));
end

X_test = X_norm(idx_test,:); y_test = y(idx_test,:);
X_cv = X_norm(idx_cv,:); y_cv = y(idx_cv,:);
X_train = X_norm(idx_train,:); y_train = y(idx_train,:);

figure; 
subplot(3,1,1); imagesc(X_test); title('Test set');
subplot(3,1,2); imagesc(X_cv); title('CV set');
subplot(3,1,3); imagesc(X_train); title('Train set');

%% 
% #####################################################################
% ==== Make the order of voxels to be introduced to the SA =====
% #####################################################################
% The normal order
featureID_sorted = [1:n];
% featureID_order = featureID_sorted(randperm(n)); % Randomly shuffle the voxel order
featureID_order = featureID_sorted ; % the MI-descend order
%% 
% #####################################################################
% #####################################################################
% ======= after this point is SA framework ===========
% #####################################################################
% #####################################################################

% #####################################################################
% Simulated Annealing optimization 
% Template file: Here we use the SA_framework_2.m
% =====================================================================
% #####################################################################


% option for display the learning curve @#$% user-defined
option = {};
option.display = 1; % (0)1 (not )display the learning curve --> 1,
option.display_final = 0; % display the final result --> 1
option.display_T = 1; % (0)1 (not )display temperature on the learning curve
option.pause = 0.01; % sec paused per each plot
option.dy_plot = 0.1; % the printing offset for temperature

% =================================================================
%                       SA framework
% =================================================================

% ===== SA operating parameters ==== @#$% user-defined
% @#$% user-defined objective function
objFunction = @ObjectiveFunction1; 
c = 0;
% @#$% user-defined classification algorithm
% SVM
PredictAccuracy = @svmPredictAccuracy; 
param = '-q -t 0 -s 0 -b 1 -c 1 -g 1';

% @#$% user-defined solution mapping function
nextSolution = @NextSolution3; 

% @#$% SA parameters
T = 1;
T_stop = 0.005; % the stopping temperature 
alpha = 0.9;
itt_max = 2000;
Rep_T_max = 100; % max number of Rep in one temperature T
Rep_accept_max =30; % if solution accepted this many, then update T
kc = 10; % k-constant. The smaller kc --> less solution accepted @#$% user-defined
% ====================================================

% ====================================================
% option for the plot/display
if option.display == 1
    figure(1919); subplot(2,1,1); hold on;
    title('learning curve'); xlabel('iterations'); ylabel('obj function');
    if isempty(option.dy_plot)
        dy_Plot = 600; % @#$% user-defined
    else
        dy_Plot = option.dy_plot; % @#$% user-defined
    end
end
% ====================================================

% ====================================================
% initialize or evaluate the solution for user-defined function  % @#$% user-defined
x_curr = zeros(1,n); x_curr(1:5) = 1; % starting with choosing first 5 voxels
F_curr = -1e+9; % F_curr = objFunction(x_curr);
x_best = x_curr;
F_best = F_curr;
storage_new = {};
% ====================================================

% ====================================================
% initial value for SA
Rep_T = 0;
Rep_accept = 0;
itt = 1;
converges = 0;
is_best_updated = 0; % toggle flag set to 1 when the best  solution is updated

% Start SA optimization
while converges == 0 && itt <= itt_max
    % pick a new solution x_new
    x_new = nextSolution(x_curr, itt, n, featureID_order); % @#$% user-defined
    
    % @#$% user-defined objective function
    % F_new = objFunction(x_new): Calculate F_new from x_new: 
    F_new = objFunction(n, c, y_cv, X_cv(:,x_new), y_train, X_train(:,x_new), param, PredictAccuracy);
    
    % Anything you want to keep is here % @#$% user-defined
    storage_new.itt = itt;
    storage_new.sol = x_new;
    storage_new.value = F_new;
    storage_new.alpha = alpha;
    storage_new.T = T;
    
    
    % Optimization
    if option.display == 1 % real-time display the learning curve
        % Check: update solution
        if F_new >= F_curr
            % back up before update
            F_prev = F_curr; % backup the previous
            x_prev = x_curr; % backup the previous
            % update
            F_curr = F_new; % accept the solution
            x_curr = x_new; % accept the solution
            Rep_accept = Rep_accept + 1;
            
            % plot the learning curve (The new solution is accepted because
            % it outperform the current solution) 
            figure(1919); subplot(2,1,1); plot(itt,F_curr,'b.','MarkerSize',20);
            
            % check the best solution || the fewer number of features
            if F_new >= F_best
%             if F_new > F_best || ( F_new==F_best && sum(x_new(:),1)<=sum(x_best(:),1) )
                F_best = F_new;
                x_best = x_new;
                storage_best = storage_new;
                % plot the learning curve (The new solution not only better
                % than the current solution, but also outperforms the best
                % solution so far)  
                figure(1919); subplot(2,1,1);  plot(itt,F_best,'go','MarkerSize',10);
                is_best_updated = 1;
                subplot(2,1,2); imagesc(x_best);
                title(['iteration: ',num2str(itt),...
                    ',   best accuracy: ',num2str(100*F_best),...
                    '%,    # voxel: ',num2str(sum(x_new(:)))]);
            end
        elseif exp((F_new-F_curr)/(kc*T)) > rand(1)
            % back up before update
            F_prev = F_curr; % backup the previous
            x_prev = x_curr; % backup the previous
            % update
            F_curr = F_new; % accept the solution
            x_curr = x_new; % accept the solution
            Rep_accept = Rep_accept + 1;
            
            % plot the learning curve (Although the new solution is poorer
            % than the current one, the new solution is accepted by "luck") 
            figure(1919); subplot(2,1,1); plot(itt,F_curr,'c.','MarkerSize',20);
        else
            % plot the learning curve (The new solution is rejected)
            figure(1919); subplot(2,1,1); plot(itt,F_curr,'rx','MarkerSize',10);
        end
        
        pause(option.pause);
        
    else % NOT displayed
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
        % plot the temperature
        if option.display == 1 && option.display_T == 1
            figure(1919); subplot(2,1,1); plot(itt,F_curr,'ks-');
            line([itt itt],[F_curr F_curr + dy_Plot],'Color','k');
            text(itt,F_curr + dy_Plot,num2str(T),'Rotation',90);
        end
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

% display the result
% fprintf('The best solution is x_best= %4.3e and F_new= %4.4e\n',x_best, F_best);
disp(['best cross-validate accuracy = ',num2str(F_best*100),'%']);
%% Best test accuracy
x_best = logical(x_best);
acc_test = objFunction(n, c, y_test, X_test(:,x_best), y_train, X_train(:,x_best), param, PredictAccuracy);
disp(['best test accuracy = ',num2str(acc_test*100),'%']);
disp(['voxels used = ',num2str(sum(x_best(:),1)/length(x_best)*100),'%']);
toc;