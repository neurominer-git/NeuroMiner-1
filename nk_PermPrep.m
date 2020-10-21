function dat = nk_PermPrep(dat, analind)

global CV CLUST

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PREPARATIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
analysis = dat.analysis{analind};
spm5ver = nk_CheckSPMver;

% ************************ LOAD PREPROCESSED DATA *************************
[preprocmat, dims] = nk_GenPreprocMaster(spm5ver, dat.id, CV);
strout = nk_Preprocess_StrCfg(CLUST, [], [], dims);

if length(dims) > 1
    fprintf('\n\nYou loaded preprocessed data with multiple dimensions.')
    fprintf('\nThe following dimensions were found:')
    fprintf('\n------------------------------------------------------')
    for i=1:length(dims)
        ind = find(analysis.dims == dims(i));
        if isempty(ind)
            fprintf('\nDim. index %g: dimensionality = %g <== NOT in selected analysis!',i,dims(i));
        else
            if dims(i)==analysis.best_dim
                bestdimind = i;
                fprintf('\nDim. index %g: dimensionality = %g => mean performances: CV1 = %1.4f, CV2 = %1.4f <== best dimensionality', ...
                    i,dims(i),analysis.TrainPerformanceMean(ind), analysis.TestPerformanceMean(ind));
            else
                fprintf('\nDim. index %g: dimensionality = %g => mean performances: CV1 = %1.4f, CV2 = %1.4f', ...
                    i,dims(i),analysis.TrainPerformanceMean(ind), analysis.TestPerformanceMean(ind));
            end
        end
    end
    dim_index = nk_input(['Select a dimension index [min = ' ...
                        num2str(1) ', max = ' ...
                        num2str(length(analysis.dims)) '] for the permutation analysis'], 0,'e',bestdimind);
else
    dim_index = 1;
end
preprocmat = preprocmat(dim_index,:,:);

% *************************** LOAD GDDATAMATs ******************************
[gdmat, gddims] = nk_GenCVdataMaster(spm5ver, dat.id, CV, dims(dim_index));

% use only those dimensions that have been defined by the user
gdmat = gdmat(gddims == dims(dim_index),:,:);

% ************* CONSTRUCT PERMUTATION ANALYSIS INPUT STRUCTURE *************
inp.label      = dat.label;             % targets to predict
if strcmp(dat.modeflag,'regression')    
    inp.nclass = 1;                     % in case of regression only 1 bin. comp.
else
    inp.nclass = length(CV.class{1,1}); % # of binary comparisons
end
inp.l           = length(dat.label);    % # of subjects
inp.cv          = CV;                   % cross-validation structure
inp.preprocmat  = preprocmat;           % paths to preprocessed data files
inp.gdmat       = gdmat;                % paths to precomputed GD structures
inp.nperms      = nk_input('# of permutations',0,'e');
inp.ovrwrt      = nk_input('Overwrite exisiting PERMdatamats?',0,'yes|no',[1,0],1);

switch analysis.params.SVM.kernel.typ
    case {' -t 0',' -t 4',' -t 5', 'lin'}
        inp.nogexp = 1;
    otherwise
        inp.nogexp = 0;
end

switch analysis.params.SVM.prog
    case 'MikRVM'
        inp.nocexp = 1;
    otherwise
        inp.nocexp = 0;
end

% ************************* CV2 GRID ACTIONS SETUP ************************
if ~exist('GridAct','var')
    [ix, jx] = size(CV.TrainInd);
    GridAct = nk_CVGridSelector(ix,jx);
end

%%%%%%%%%%%%%%%%%%%%%%% RUN PERMUTATION ANALYSIS %%%%%%%%%%%%%%%%%%%%%%%%%%% 
for i=1:length(dim_index) % Loop through selected dimensions
    inp.dim_index = i;
    inp.dimension = dims(dim_index(i));        % current dimension to analyze
    ind = find(analysis.dims ==  dims(dim_index(i)));
    if isempty(ind), fprintf('\nNo GDdims structure found for selected analysis.'), continue; end
    inp.analysis = analysis.GDdims{ind};
    stranalysis = ['_dim' num2str(inp.dimension) '_p' num2str(inp.nperms)];
    analysis.Perm{i} = nk_PermX(inp, stranalysis, dat.id, GridAct);
end

dat.analysis{analind} = analysis;

return
