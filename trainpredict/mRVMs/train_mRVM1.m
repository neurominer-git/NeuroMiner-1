%% mRVM1
% authors: Ioannis Psorakis psorakis@gmail.com, Theodoros Damoulas theo@dcs.gla.ac.uk
%    Copyright NCR, 2009
%
%    Use of this code is subject to the Terms and Conditions of the Research License Agreement,
%    agreed to when this code was downloaded, and a copy of which is available at
%    http://www.dcs.gla.ac.uk/inference/pMKL/Terms_and_Conditions.html.
% INPUTS:
% =======
% USE '-i' as input_flag if you want to provide the inputs interactively
% through the console
% OR
% use '-p' as input_flag if you want to provide the inputs as function
% parameters (see function declaration above)
% In case you use '-p' please provide the following:
% X: of size N x D (N samples, D features) is the training set.
% t: of size C x N (C classes, N samples) is the training labels.
% standardize_flag: [boolean] turns data standardization ON/OFF
% convergence_used: values [1 or 2] is the training termination criterion (see
%                   conv.1 and conv.2 of theoretical background
% kernel_type: string can be either 'gaussian', or 'polynomial' or 'linear;
% kernel_param: for linear kernel put any value
% plot_flag: [1 or 0] plots the number of relevant vectors during training
% dataset_name: auxiliary label
% In case you use '-i' please note:
% each dataset file must contain the necessary variables, which much be in the correct format. Those are:
% 
% 1 the class labels, say "t". This variable must be a C X N dimensional array
%   where C is the number of classes and N is the size of the data. Say we have a
%   problem with 100 samples and 3 classes. The 5th sample would belong to the 2nd
%   class if t(:,5)' = [0 1 0].
%   For datasets which do have independent training and test sets, there should be two of these
%   variables. E.g. tTrain and tTest
% 
% 2 the data, say "X". This variable should be a N X D where N is the number of samples
%   and D is the number of features. For datasets which do have independent training
%   and test sets, there should be two of there variables. E.g. Xtrain Xtest. Also, for multi-kernel problems
%   there should be one such variable for each feature space.
% OUTPUTS:
% =======
% OUTPUT is an object that has a range of properties:
% model_used: the name of the algorithm (e.g mRVM1);
% dataset_name: the name of the dataset;
% N_total: the total number of samples in the original dataset;
% N_prototypical: the number of relevance vectors extracted from the algorithm;
% 
% X_prototypical: the relevance vectors (or prototypical samples) extracted from the algorithm;
% 
% X_prototypical_standardized: same as above, but standardized;
% K_sparse: the sparse training kernel, that uses only the relevance vectors;
% W_sparse: the sparse regressors matrix, that uses only the relevance vectors;
% 
% active_sample_original_indices: the original indices of the relevance vectors in the dataset;
% 
% sources: number of kernels;
% b: kernel mixing coefficients;
% 
% kernel_type: kernel(s) used;
% kernel_param: kernel parameter(s) used;

function OUTPUT = train_mRVM1(input_flag,X,t,standardize_flag,convergence_used,Nmax,kernel_type,kernel_param,plot_flag,dataset_name)

global VERBOSE
%disp('---------------- mRVM1 ----------------')
model_used = 'mRVM1';

if ~exist('input_flag','var')
    input_flag = '-i';
end

if strcmp(input_flag,'-i')
    %% INTERFACE
    % provide the dataset file name
    dataset_name = input('please enter the file name of the dataset > ','s');
    
    % attempt to open the file name, if not, stop
    try
        load(dataset_name);
    catch ME
        ME.stack;
        error('No dataset exists with that name, exiting...')
    end
    
    X = input('please provide the training set variable name > ');
    
    standardize_flag = logical(input('standardize data? (1/0) > '));
    
    % KERNEL SETUP
    
    fprintf(1,'Please enter the kernel number:\n1. Gaussian\n2. Polynomial\n3. Linear\n');
    kernel_ID = input('your choice > ');
    switch kernel_ID
        case 1
            kernel_type='gaussian';
            kernel_param = input('please enter the kernel parameter > ');
        case 2
            kernel_type = 'polynomial';
            kernel_param = input('please enter the kernel parameter > ');
        case 3
            kernel_type = 'linear';
            kernel_param='';
    end
    
    t = input('please provide the train labels variable name > ');
    
    %%%
    fprintf('Please select convergence criterion:\n 1. trivial change in hyperpriors \n 2. trivial change in hyperpriors and minimum number of iterations N \n 3. Maximum iterations\n');
    convergence_used = input('choice (1 or 2 or empty if not sure) > ');
    
    if isempty(convergence_used)
        convergence_used = 2;
        Nmax = size(X,1);
    elseif convergence_used==1
        Nmax = -1;
    elseif (convergence_used == 2) || (convergence_used == 3)
        Nmax = input('Enter maximum number of iterations (leave empty for Ntrain) > ');
        if isempty(Nmax)
            Nmax = size(X,1);
        end
    end
    
    plot_flag = input('plot relevant vectors progresion? (1/0) > ');
    

%elseif strcmp(input_flag,'-f')    
%elseif strcmp(input_flag,'-p')    
end

%% INITIALIZATIONS
if standardize_flag
    Xtrain = standardize_data(X);
else
    Xtrain = X;
end

N = size(Xtrain,1);
D = size(Xtrain,2);
C = size(t,1);

% initialize the auxiliary variables Y to follow the target labels of the
% training set
Y = 10*rand(C,N).*t + rand(C,N);
% set all scales to infinity
A = inf*ones(N,C);
% set all regressor values to zero
W = zeros(N,C);

if plot_flag
    number_of_RVs = zeros(N+1,1);
end

if isempty(Nmax)
    Nmax = size(X,1);
end

% mRVM1 starts with an empty model
active_samples=[];
%% PERFORM KERNEL SETUP
%disp('building kernel(s)...')
Ktrain = build_standarized_kernel...
    (Xtrain,Xtrain,kernel_type,kernel_param);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SELECT FIRST SAMPLE
% here we follow an informative methodology for the selection of the first
% sample to include in the model. We select the sample i, for which k_i has
% the largest projection to the auxiliary variables Y normalized by
% ||k_||^2 and the number of classes C
Initial_CON = sqrt(sum((Ktrain*Y').^2,2)) ./ (C*sum(Ktrain.^2,2));
current_sample_index = find(Initial_CON == max(Initial_CON),1);

% add this sample to the model.
active_samples = [active_samples;current_sample_index];

%% BUILD KERNEL, UPDATE K*,ki,A
% take the i-th row of the training Kernel for that sample
ki = Ktrain(current_sample_index,:);
% If at any point we have M samples in the model, Kstar is MXN. In our case
% where M=1, Kstar=ki
Kstar = Ktrain(active_samples,:);

% update the scales for the selected sample.
A(current_sample_index,:) = update_a_mRVM1_first_time(ki,Y);
% Astar is MXN for M samples in the model. Initially Astar = 1XN because we
% have only one sample in the model.
Astar = diag(A(active_samples,1));
%% LOOP
% Estimate the quantity KKA_inv which will be used for the calculation of
% the contribution of a sample to the model. It shall be also used for the
% update of the W posteriors.
KKA_inv = chol_inv(Kstar*Kstar' + Astar);

% initiallize some iterators and auxiliary variables
mRVM1_iterator=0;
La=1;

% set all convergence flags to false.
converged_loga_N=false;
converged_loga=false;

%disp('Training...')
tic
while true
    %% Results of previous iteration
    %   store some results from the previous iteration.
    Aold=diag(Astar);
    Laold=La;
    Wold = W;
    Yold = Y;
    %% Ki
    % calculate the kernel function for the selected sample.
    ki = Ktrain(current_sample_index,:);
    %% CHECK IF ALREADY IN
    % check if the selected sample is already inside the relevant vectors
    % collection
    is_already_in = ...
        logical(sum(1*ismember(active_samples,current_sample_index)));
    % get its scale
    ai = A(current_sample_index,:);
    %% CALCULATE Sm,Qcm
    % calculate the sparsity and quality factors (See theoretical
    % background)
    if mRVM1_iterator==0
        Sm = ki*ki' - ki*Kstar'*KKA_inv*Kstar*ki';
        Qcm = ki*Y' - ki*Kstar'*KKA_inv*Kstar*Y';
        %% UPDATE s,q based on sample existence
        if ~is_already_in
            si = Sm;
            qci = Qcm;
        else
            si = (ai(1)*Sm)/(ai(1)-Sm);
            qci = (ai(1)*Qcm)./(ai(1)-Sm);
        end
        %% calculate sum of qci^2
        qsum2=sum(qci.^2);
        
        does_contribute = qsum2 > C*si;
    else
        Sm = BIG_SM(current_sample_index);
        Qcm = BIG_QM(current_sample_index,:);
        
        if ~is_already_in
            si = Sm;
            qci = Qcm;
        else
            si = (ai(1)*Sm)/(ai(1)-Sm);
            qci = (ai(1)*Qcm)./(ai(1)-Sm);
        end
        qsum2=sum(qci.^2);
        does_contribute = qsum2 > C*si;
    end
    %% CHECK CONDITIONS
    % based on the contribution value of the current sample and its
    % existence or not in the relevant vector collection, for each case we
    % do the following:
    if does_contribute && is_already_in
        % if the sample does contribute to the model solution and it is
        % already in the model keep it and update its scales a
        A(current_sample_index,:) = update_a_mRVM1(C,si,qci);
        
        % update the matrices
        Astar = diag(A(active_samples,1));
        KKA_inv = chol_inv(Kstar*Kstar' + Astar);
        
        % there has been a change in the model (due to the scale update so
        % we need to update the posteriors).
        update_posterior_flag=true;
        
    elseif does_contribute && ~is_already_in
        % is the sample can contribute to the model, but it is not inside
        % it, we include it.
        A(current_sample_index,:) = update_a_mRVM1(C,si,qci);
        
        active_samples = [active_samples current_sample_index];
        
        active_samples = sort(active_samples);
        
        % update matrices
        Astar = diag(A(active_samples,1));
        Kstar = Ktrain(active_samples,:);
        KKA_inv = chol_inv(Kstar*Kstar' + Astar);
        
        % as the model changed due to the inclusion of a sample, we need to
        % update the posteriors
        update_posterior_flag=true;
        
    elseif ~does_contribute && is_already_in && length(active_samples)~=1
        % if the sample does not contribute to the model solution and it is
        % in the model, we remove it
        A(current_sample_index,:) = inf;
        W(current_sample_index,:)=0;
        active_samples(active_samples==current_sample_index)=[];
        
        active_samples = sort(active_samples);
        
        % update matrices
        Astar = diag(A(active_samples,1));
        Kstar = Ktrain(active_samples,:);
        KKA_inv = chol_inv(Kstar*Kstar' + Astar);
        
        % as we removed one sample in the model, we have to update our
        % posteriors.
        update_posterior_flag=true;
    else
        % the posterior update is not needed if the current sample does not
        % contribute to the solution and it is not in the model.
        update_posterior_flag=false;
    end
    non_active=1:1:N;
    non_active(active_samples)=[];
    %% UPDATE POSTERIORS
    % if there was a change in the model, update the regressors W and the
    % auxiliary variables Y.
    if update_posterior_flag
        %% W
        W(active_samples,:) = KKA_inv * Kstar * Y';
        %% Y
        F=Kstar'*W(active_samples,:);
        for c=1:C
            pos=find(t(c,:)==1);
            [Z, Yt(pos,:)] = YTruncate(F(pos',:),c,1000);
        end
        Y=Yt';
        %% BIG MATRICES
        % those matrices are used to evaluate the contribution of each
        % sample of the training set. We initialize:
        %non_active=1:1:N;
        %non_active(active_samples)=[];
        BIG_SI = zeros(1,N);
        BIG_QCI = zeros(N,C);
        
        % BIG_SI stores the `sparsity factor' (see theoretical background)
        % for each sample in the training set. BIG_QCI stores the `quality
        % factor' (see theoretical background) for each sample in the model
        % and across classes.
        
        % BIG_M and BIG_QM are auxiliary tables for the calculation of the
        % above factors a different update scheme is followed for included
        % and not-included samples of the model (see theoretical
        % background)
        BIG_SM = diag(Ktrain*Ktrain - Ktrain*Kstar'*KKA_inv*Kstar*Ktrain); % 1XN
        BIG_QM = (Ktrain*Y' - Ktrain*Kstar'*KKA_inv*Kstar*Y'); % NXC
        
        BIG_SI(non_active) = BIG_SM(non_active);
        BIG_QCI(non_active,:) = BIG_QM(non_active,:);
        
        BIG_SI(active_samples) =...
            (A(active_samples,1).*BIG_SM(active_samples))./...
            (A(active_samples,1)-BIG_SM(active_samples));
        
        BIG_QCI(active_samples,:) =...
            (A(active_samples,:).*BIG_QM(active_samples,:))./...
            repmat(A(active_samples,1)-BIG_SM(active_samples),1,C);
        
        % having the sparsity and quality factor for each individual
        % sample, we calculate each sample contribution, which is stored in
        % to the BIG_CON matrix. This matrix shall be used for the
        % informative sample selection scheme (see theoretical background).
        BIG_CON = sum(BIG_QCI.^2,2) - C*BIG_SI';
    end
    %% UPDATE CONVERGENCE
    % Calculations for the estimation of the convergence based on the
    % change of the scales matrix A.
    Acurr=diag(Astar);
    if length(Acurr)>length(Aold)
        Aold_aux = [Aold;eps*ones(length(Acurr)-length(Aold),1)];
        Acurr_aux=Acurr;
    elseif length(Aold)>length(Acurr)
        Acurr_aux =[Acurr;eps*ones(length(Aold)-length(Acurr),1)];
        Aold_aux=Aold;
    else
        Aold_aux = Aold;
        Acurr_aux = Acurr;
    end
    
    % LOGA convergence
    loga = abs(log(Acurr_aux)-log(Aold_aux));
    
    
    %Convergence before N
    are_loga_all_small = prod(1*(loga<1e-2));
    are_non_active_all_negative = prod(1*(BIG_CON(non_active)<=0));
    
    if ~converged_loga
        converged_loga = are_loga_all_small &&...
            are_non_active_all_negative;
    end
    
    %Convergence after N
    if ~converged_loga_N
        converged_loga_N = are_loga_all_small &&...
            are_non_active_all_negative && mRVM1_iterator>Nmax;
    end
    
    % stop when converged
    if (convergence_used == 1 && converged_loga) || (convergence_used == 2 && converged_loga_N) ||...
            (convergence_used == 3 && mRVM1_iterator == Nmax)
        break;
    end
    
    % break if not converged after 10*N
    if mRVM1_iterator == 10 * N
        break;
        warning('not converged - training terminated after maximum iterations')
    end
    %% PLOTS
    if plot_flag
        number_of_RVs(mRVM1_iterator+1) = length(active_samples);
        plot(number_of_RVs(1:mRVM1_iterator+1));
        title('Relevant vectors')
        drawnow;
    end
    %% INFORMATIVE SAMPLE SELECTION
    % based on the contributions matrix BIG_CON, perform an informative
    % selection scheme for the next sample which will be evaluated. At
    % first, we try to find in a non-included sample i can contribute to
    % the model solution (BIG_CON(i)>0). If more than one sample exists
    % with that criterion, select the one with the maximum contribution. If
    % no such sample exists, for each sample inside the model check if
    % there is one which has a negative contribution. If more than one such
    % samples exist, select the one with the biggest by absolute value
    % negative contribution.
    
    % perform some intermediate calculations
    BIG_CONtmp = BIG_CON;
    BIG_CONtmp(active_samples)=-inf;
    
    non_active=1:1:N;
    non_active(active_samples)=[];
    BIG_CONtmp2 = BIG_CON;
    BIG_CONtmp2(non_active)=inf;
    
    % and apply the criteria described above
    if ~isempty(BIG_CONtmp(BIG_CONtmp>0))
        % ARE THERE ANY NON-ACTIVE WHICH CAN CONTRIBUTE TO THE MODEL?
        current_sample_index = ...
            find(BIG_CONtmp==max(BIG_CONtmp(BIG_CONtmp>0)),1);
        %disp('Inactive which can contribute selected:')
        %current_sample_index
    elseif sum(1*(BIG_CON(active_samples)<=0))~=0
        % ARE THERE ANY ACTIVE WHICH DO NOT CONTRIBUTE?
        current_sample_index = ...
            find(BIG_CONtmp2==min(BIG_CONtmp2),1);
        %disp('Active which did not contribute selected')
        %current_sample_index
    else
        % UPDATE THE A OF ONE RANDOM ACTIVE SAMPLE
        current_sample_index = ...
            active_samples(ceil(1+rand*(length(active_samples)-1)));
        %disp('Random active sample selected') current_sample_index
    end
    %% INCREMENT
    mRVM1_iterator=mRVM1_iterator+1;
end
%disp('TRAINING FINALIZED.');
%disp('-------------------');
%toc
tElapsed = toc;
if ~VERBOSE, fprintf('\tMKL-RVM1: completed in %1.2f seconds (%1.1f%% active samples).',tElapsed, length(active_samples)*100/N); end
%fprintf('Total iterations until convergence: %d\n',mRVM1_iterator)
%fprintf('Number of prototypical samples: %d out of %d total.\n',length(active_samples),N)

if ~exist('dataset_name','var')
    dataset_name = 'untitled_run';
end

% 
if standardize_flag
    Xtrain = Xtrain(active_samples,:);
else
    Xtrain = nan;
end

OUTPUT = mRVM_train_output(model_used, dataset_name, N, X(active_samples,:),standardize_flag, Xtrain, Ktrain(active_samples,active_samples),...
    W(active_samples,:), active_samples,kernel_type,kernel_param);
end

%% AUXILIARY FUNCTIONS
% authors: Ioannis Psorakis psorakis@gmail.com, Theodoros Damoulas
% theo@dcs.gla.ac.uk
%    Copyright NCR, 2009 
%
%    Use of this code is subject to the Terms and Conditions of the Research License Agreement,
%    agreed to when this code was downloaded, and a copy of which is available at 
%    http://www.dcs.gla.ac.uk/inference/pMKL/Terms_and_Conditions.html.

function ai = update_a_mRVM1(C,si,qci)

result = (C*si^2) / (sum(qci.^2)-C*si);

if result<1e-5
    result = 1e-6;
end

ai=result;
end

function a_new = update_a_mRVM1_first_time(ki,Y)

C = size(Y,1);

nominator = sum(ki.^2);

kis = repmat(ki,C,1);

denominator = sum(sum(kis.*Y,2).^2) / (C*nominator -1);

a_new = nominator/denominator;
end