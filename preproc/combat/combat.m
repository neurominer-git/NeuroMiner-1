function [bayesdata,estimators] = combat(dat, batch, mod, estimators)

% This function was originally sent to Dom Dwyer from Russel Shinohara and
% is based on the following papers and the git repo below:
% Jean-Philippe Fortin, Drew Parker, Birkan Tunc, Takanori Watanabe, Mark A Elliott, Kosha Ruparel, David R Roalf, Theodore D Satterthwaite, Ruben C Gur, Raquel E Gur, Robert T Schultz, Ragini Verma, Russell T Shinohara. Harmonization Of Multi-Site Diffusion Tensor Imaging Data. NeuroImage, 161, 149-170, 2017
% Jean-Philippe Fortin, Nicholas Cullen, Yvette I. Sheline, Warren D. Taylor, Irem Aselcioglu, Philip A. Cook, Phil Adams, Crystal Cooper, Maurizio Fava, Patrick J. McGrath, Melvin McInnis, Mary L. Phillips, Madhukar H. Trivedi, Myrna M. Weissman, Russell T. Shinohara. Harmonization of cortical thickness measurements across scanners and sites. NeuroImage, 167, 104-120, 2018
% https://github.com/Jfortin1/ComBatHarmonization
%
% INPUTS
% dat        - (required) a data matrix (p x n) for which the p rows are features, and the n columns are participants. 
% batch      - (required) a numeric or character vector of length n indicating the site/scanner/study id.
% mod        - a (n x p) matrix containing the outcome of interest (e.g., illness) and other covariates that will not be removed (i.e., they will be statistically preserved). 
% estimators - stored parameters from a previous combat run used to apply to test data
%
% OUTPUTS
% bayesdata  - the corrected data matrix
% estimators - parameters used to correct data that can be applied to new data
%
% For directions on function use in the training case, see the README.md
% file that was produced by the developers. 
%
% The modification was to allow external validation on test data. I did
% this by creating an "estimators" structure that is output when the
% algorithm is run on the training data. This can then be entered as a
% field when the test data is run. 
%
% When running on training data, do not include the "estimators"
% variable and run with: 
% [bayesdata,estimators] = combat(dat, batch, mod)
%
% When running on test data, include the "estimators" field that comes from
% the training data: 
% [bayesdata] = combat(dat, batch, mod, estimators)
%
% This is the same logic as regression, so it means that if a test site is
% not included in the training sites then it won't work. 
%
% Original function: JP Fortin?, don't know when it was produced. Need to follow-up.  
% Edits: Dom edits to add external validation 1July2020.  
%%

switch nargin
    
    case 1 
        fprintf('\nDefine batch variable')
        return
        
    case 2
        fprintf('\nDefine mod as [] if no mod exists')
        return
        
    case 3
        
        % Check data and define batches
        [sds] = std(dat')';
        wh = find(sds==0);
        [ns,ms] = size(wh);
        if ns>0
            error('Error. There are rows with constant values across samples. Remove these rows and rerun ComBat.')
        end
        batchmod = dummyvar(batch);
        Ix = sum(batchmod)==0; batchmod(:,Ix)=[]; 
        n_batch = size(batchmod,2);     
        levels = unique(batch);         
        fprintf('[combat] Found %d batches\n', n_batch);
        
        batches = cell(0);
        for i=1:n_batch
            batches{i}=find(batch == levels(i));
        end
        n_batches = cellfun(@length,batches);
        n_array = sum(n_batches);
        
        % Creating design matrix and removing intercept:
        design = [batchmod mod];
        intercept = ones(1,n_array)';
        wh = cellfun(@(x) isequal(x,intercept),num2cell(design,1));
        bad = find(wh==1);
        design(:,bad)=[];
        
        
        fprintf('[combat] Adjusting for %d covariate(s) of covariate level(s)\n',size(design,2)-size(batchmod,2))
        % Check if the design is confounded
        if rank(design)<size(design,2)
            nn = size(design,2);
            if nn==(n_batch+1)
                error('Error. The covariate is confounded with batch. Remove the covariate and rerun ComBat.')
            end
            if nn>(n_batch+1)
                temp = design(:,(n_batch+1):nn);
                if rank(temp) < size(temp,2)
                    error('Error. The covariates are confounded. Please remove one or more of the covariates so the design is not confounded.')
                else
                    error('Error. At least one covariate is confounded with batch. Please remove confounded covariates and rerun ComBat.')
                end
            end
        end
        
        
        fprintf('[combat] Standardizing Data across features\n')
        B_hat = inv(design'*design)*design'*dat';
        
        %Standarization Model
        grand_mean = (n_batches/n_array)*B_hat(1:n_batch,:);
        var_pooled = ((dat-(design*B_hat)').^2)*repmat(1/n_array,n_array,1);
        stand_mean = grand_mean'*repmat(1,1,n_array);
        
        if not(isempty(design))
            tmp = design;
            tmp(:,1:n_batch) = 0;
            stand_mean = stand_mean+(tmp*B_hat)';
        end
        s_data = (dat-stand_mean)./(sqrt(var_pooled)*repmat(1,1,n_array));
        
        %Get regression batch effect parameters
        fprintf('[combat] Fitting L/S model and finding priors\n')
        batch_design = design(:,1:n_batch);
        gamma_hat = inv(batch_design'*batch_design)*batch_design'*s_data';
        delta_hat = [];
        for i=1:n_batch
            indices = batches{i};
            delta_hat = [delta_hat; var(s_data(:,indices)')];
        end
        
        %Find parametric priors:
        gamma_bar = mean(gamma_hat');
        t2 = var(gamma_hat');
        delta_hat_cell = num2cell(delta_hat,2);
        a_prior=[]; b_prior=[];
        for i=1:n_batch
            a_prior=[a_prior aprior(delta_hat_cell{i})];
            b_prior=[b_prior bprior(delta_hat_cell{i})];
        end
        
        fprintf('[combat] Finding parametric adjustments\n')
        
        gamma_star =[]; delta_star=[];
        for i=1:n_batch
            indices = batches{i};
            temp = itSol(s_data(:,indices),gamma_hat(i,:),delta_hat(i,:),gamma_bar(i),t2(i),a_prior(i),b_prior(i), 0.001);
            gamma_star = [gamma_star; temp(1,:)];
            delta_star = [delta_star; temp(2,:)];
        end
        
        fprintf('[combat] Adjusting the Data\n')
        bayesdata = s_data;
        j = 1;
        for i=1:n_batch
            indices = batches{i};
            bayesdata(:,indices) = (bayesdata(:,indices)-(batch_design(indices,:)*gamma_star)')./(sqrt(delta_star(j,:))'*repmat(1,1,n_batches(i)));
            j = j+1;
        end
        bayesdata = (bayesdata.*(sqrt(var_pooled)*repmat(1,1,n_array)))+stand_mean;
        
        % save out estimators for external application
        estimators.levels = levels;
        estimators.n_batch = n_batch;
        estimators.var_pooled = var_pooled;
        estimators.B_hat = B_hat;
        estimators.grand_mean = grand_mean;
        estimators.gamma_star = gamma_star;
        estimators.delta_star = delta_star;
        
    case 4
        
        if isstruct(estimators) && ~isempty(estimators)
            [bayesdata] = combat_test(dat, batch, mod, estimators);
        else
            fprintf('Estimators is empty. Enter estimators structure or do not enter argument')
        end
        
end

end


function [bayesdata] = combat_test(dat, batch, mod, estimators)

% This function takes the estimators from a training analysis and applies
% them to the test data. In the script above, the data is first
% standardized using a weighted grand mean of the data divided by a measure
% of variance. Coefficients are then generated to ultimately get the
% gamma_star and delta_star matrices that are used to correct the
% standardized data. The logic of this script was to use the grand mean,
% weights for the grand mean, variance, and the gamma/delta star from the
% training analysis and then apply them to the test data. 

% separate section for test data
levels = estimators.levels;
n_batch = estimators.n_batch;
var_pooled = estimators.var_pooled;
B_hat = estimators.B_hat;
grand_mean = estimators.grand_mean;
gamma_star = estimators.gamma_star;
delta_star = estimators.delta_star;

% Initialize the batchmod and populate
n_array = size(dat,2);
batchmod = zeros(n_array,n_batch);

% Can't use dummyvar to populate when only first group is present so loop
% to populate the batchmod
for i=1:n_batch
    batchmod(:,i) = batch==levels(i); 
end

fprintf('[combat] Found %d batches\n', n_batch); % sanity check

batches = cell(0);
for i=1:n_batch
    batches{i}=find(batch == levels(i)); % get the indices for each subject. Some batches{i} will be empty
end

n_batches = cellfun(@length,batches);

% Creating design matrix and removing intercept:
design = [batchmod mod];

% This section is from the script above. 
% In test circumstances I can't remove columns that are all a single number because of the
% case where there is only one site. This step was probably to reduce
% redundancy--in their words "counfounded matrices". 
% intercept = ones(1,n_array)';
% wh = cellfun(@(x) isequal(x,intercept),num2cell(design,1));
% bad = find(wh==1);
% design(:,bad)=[];

batch_design = design(:,1:n_batch);

fprintf('[combat] Standardizing the Data Based on Training Sample\n')

stand_mean = grand_mean'*repmat(1,1,n_array);

if not(isempty(design))
    tmp = design;
    tmp(:,1:n_batch) = 0;
    stand_mean = stand_mean+(tmp*B_hat)';
end

s_data = (dat-stand_mean)./(sqrt(var_pooled)*repmat(1,1,n_array));

fprintf('[combat] Adjusting the Data\n')

bayesdata = s_data;
j = 1;
for i=1:n_batch
    indices = batches{i};
    bayesdata(:,indices) = (bayesdata(:,indices)-(batch_design(indices,:)*gamma_star)')./(sqrt(delta_star(j,:))'*repmat(1,1,n_batches(i)));
    j = j+1;
end
bayesdata = (bayesdata.*(sqrt(var_pooled)*repmat(1,1,n_array)))+stand_mean;

end
