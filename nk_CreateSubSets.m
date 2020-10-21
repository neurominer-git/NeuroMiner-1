function SubSets = nk_CreateSubSets(Y)
global CV RFE MULTI MODEFL MULTILABEL VERBOSE PREPROC 

if VERBOSE && RFE.Filter.flag, fprintf('\n\nCreate feature subsets'); end

[nperms, nfolds, nvar] = size(Y.Tr);
nl = MULTILABEL.dim;

SubSets = cell(nl,1);
W = cell(nl,1);
if isfield(CV(1),'class')
    nc = numel(CV(1).class{1,1});
else
    nc = 1;
end
if MULTI.flag && MULTI.train 
    % Preprocessing has been performed in multi-group mode, otherwise
    % loop through dichotomizers in multi-group mode
    if ~RFE.Filter.binmode && ~iscell(Y.Tr{1,1,1}), 
        nc = 1; 
    end    
end
if iscell(PREPROC)
    BINMOD = PREPROC{1}.BINMOD;
else
    BINMOD = PREPROC.BINMOD;
end

for curlabel=1:nl % Label loop
    
    SubSets{curlabel} = cell(nperms, nfolds, nc, nvar);
    W{curlabel} = cell(nperms, nfolds, nc, nvar);
    
    for v=1:nvar % Variate loop
        for curclass=1:nc % Binary predictor loop
            Mx = [];
            for i=1:nperms % Inner permutation loop
                for j=1:nfolds % Inner fold loop                
                    switch RFE.Filter.binmode 
                        % Multi-group Feature ranking
                        case 0
                            if VERBOSE && RFE.Filter.flag
                                switch MODEFL
                                    case 'classification'
                                        if nc > 1
                                            fprintf('\nFilter Label %g, Variate %g => CV1 [%g, %g, Multi-Group (%s)]:', ...
                                                curlabel, v,i,j,CV.class{1,1}{curclass}.groupdesc)
                                        else
                                            fprintf('\nFilter Label %g, Variate %g => CV1 [%g, %g, Multi-Group]:',curlabel,v,i,j)
                                        end
                                    case 'regression'
                                        fprintf('\nFilter Label %g, Variate %g => CV1 [%g, %g, Regression]:',curlabel,v,i,j)
                                end
                            end
                            % Get Data
                            if ~BINMOD
                                Tr = Y.Tr{i,j,v}; 
                                Cv = Y.CV{i,j,v};
                            else
                                Tr = Y.Tr{i,j,v}{curclass}; 
                                Cv = Y.CV{i,j,v}{curclass}; 
                            end
                            % Get label
                            if MULTI.flag
                                TrL = Y.mTrL{i,j}; CvL = Y.mCVL{i,j}(:,curlabel);
                            else
                                TrL = Y.TrL{i,j}{1}; CvL = Y.CVL{i,j}{1}(:,curlabel);
                            end
                        % Binary Feature ranking
                        case 1
                            if VERBOSE && RFE.Filter.flag
                                switch MODEFL
                                    case 'classification'
                                        fprintf('\nFilter Label %g, Modality %g => CV1 [%g, %g, %s]:',curlabel, v,i,j,CV.class{1,1}{curclass}.groupdesc)
                                    case 'regression'
                                        fprintf('\nFilter Label %g, Modality %g => CV1 [%g, %g, Regression]:',curlabel, v,i,j)
                                end
                            end
                            % Get Data
                            if ~BINMOD
                                Tr = Y.Tr{i,j,v}(Y.TrInd{i,j}{curclass},:); 
                                Cv = Y.CV{i,j,v}(Y.CVInd{i,j}{curclass},:);
                            else
                                Tr = Y.Tr{i,j,v}{curclass}(Y.TrInd{i,j}{curclass},:); 
                                try
                                    Cv = Y.CV{i,j,v}{curclass}(Y.CVInd{i,j}{curclass},:); 
                                catch
                                    fprintf('problem')
                                end
                            end
                            % Get label
                            TrL = Y.TrL{i,j}{curclass}; CvL = Y.CVL{i,j}{curclass}(:,curlabel);
                    end
                    [SubSets{curlabel}{i,j,curclass, v}, ...
                        W{curlabel}{i,j,curclass,v}, Mx] = CreateSubSet(Tr, TrL(:,curlabel), Cv, CvL, Mx, [], curclass);
                end
            end
            % Computer optimum parameters across all CV1 data partitions
            % (similar to Feature generation step) for IMRelief algorithm
            if RFE.Filter.type == 6 && ~isempty(Mx)

                if numel(RFE.Filter.imrelief) > 1
                    sigma = RFE.Filter.imrelief{curclass}.sigma;
                    lambda = RFE.Filter.imrelief{curclass}.lambda;
                    nsigma = numel(sigma);       
                else
                    sigma = RFE.Filter.imrelief{1}.sigma;
                    lambda = RFE.Filter.imrelief{1}.lambda;
                    nsigma = numel(sigma);
                end

                Mx = Mx./(nperms*nfolds);

                if VERBOSE, fprintf('\nComputed CV grid of sigma & lambda across CV1 partitions:'); end
                for i=1:nsigma,
                    fprintf('\n')
                    fprintf('\t%1.1f',Mx(i,:)); 
                end

                mxcv = max(Mx(:));
                [xpos, ypos] = find(Mx == mxcv);
                if VERBOSE, fprintf('\nSelected sigma = %1.1f, lambda = %1.3f, Perf = %1.2f', ...
                    sigma(xpos(1)), ... 
                    lambda(ypos(1)), mxcv); 
                end

                % Compute average weight vector over CV1 partitions
                weights = zeros(size(Tr,2),1);
                for i=1:nperms % Inner permutation loop
                    for j=1:nfolds % Inner fold loop 
                        weights = weights + W{curlabel}{i,j,curclass,v}(:,xpos(1),ypos(1))./(nperms*nfolds);
                    end
                end
                [~, ind] = sort(weights,'descend');

                % Recompute structured feature subspaces based on average
                % weight vector 
                for i=1:nperms % Inner permutation loop
                    for j=1:nfolds % Inner fold loop 
                        if BINMOD
                            Tr = Y.Tr{i,j,v}{curclass}(Y.TrInd{i,j}{curclass},:); 
                            Cv = Y.CV{i,j,v}{curclass}(Y.CVInd{i,j}{curclass},:); 
                        else
                            Tr = Y.Tr{i,j,v}(Y.TrInd{i,j}{curclass},:); 
                            Cv = Y.CV{i,j,v}(Y.CVInd{i,j}{curclass},:);
                        end
                        % Get Label
                        TrL = Y.TrL{i,j}{curclass}(:,curlabel); CvL = Y.CVL{i,j}{curclass}(:,curlabel);
                        % Get Subspace evaluation
                        SubSets{curlabel}{i,j,curclass, v} = CreateSubSet(Tr, TrL, Cv, CvL, [], ind, curclass);
                    end    
                end
            end
        end
    end
end

% _________________________________
function [S, weights, aMx] = CreateSubSet(Y, label, Ynew, labelnew, aMx, ind, curclass)

global RFE VERBOSE

kFea 	= size(Y,2); weights = [];
if ~exist('aMx','var'), aMx = []; end

if ~RFE.Filter.flag
    
    % No ranking, No subspace learning
    S=(1:kFea)';
    
elseif ~RFE.Filter.SubSpaceFlag
    
    % Features are ranked and thresholded at a user-defined threshold
    S=(1:kFea)';
    if ~exist('ind','var') || isempty(ind)
        [ weights, ind ] = apply_filter(Y,label, kFea, Ynew, labelnew, aMx, curclass);
    end

    if RFE.Filter.RankThresh > 1
        RankThresh = RFE.Filter.RankThresh/100;
    else
        RankThresh = RFE.Filter.RankThresh;
    end
    selnum = kFea - ceil(kFea*RankThresh);
    S = S(ind(selnum));

    
    if VERBOSE, fprintf(' ranking.'); end
    
else
    % Ranked feature subspaces are created, possibly giving rise to
    % classifier ensembles (alternatively, probabilistic feature extraction
    % might be used)
    k = kFea;
    km = k-(RFE.Filter.MinNum-1);
    
    if ~exist('ind','var') || isempty(ind)
     [weights, ind, S ] = apply_filter(Y,label, k, Ynew, labelnew, aMx, curclass);
    end
    
    if RFE.Filter.type ~= 11

        S = zeros(k,km);

        while km >= 1 % Generate a number of feature subspaces
            switch RFE.Filter.type
                case 0
                    kind =  1:kFea; % Include all features
                    sx = kind;
                case {1, 2, 3, 4, 5, 6, 7, 10, 12, 13, 14}
                    kind = ind(1:k);
                    sx = 1:k;  
                case 9 % Generate unfiltered feature subspaces of increasing cardinality
                    kind = 1:k;
                    sx = 1:k;
                case 11
                    kind = 1:kFea;
                    sx = ind(:,km);    
            end
            S(kind, km) = sx;
            k=k-1;
            km=km-1;
        end
    end
end

%Convert to uintX to save space
if max(S(:)) > intmax('uint32')
    S = uint64(S);
elseif max(S(:)) > intmax('uint16')
    S = uint32(S);
elseif max(S(:)) > intmax('uint8')
    S = uint16(S);
else
    S = uint8(S);
end

function [weights, ind, S] = apply_filter(Y,label, k, Ynew, labelnew, aMx, curclass)
global VERBOSE RFE

weights=[]; ind = []; S=[]; 

switch RFE.Filter.type

    case 1 % FEAST

        ind = nk_FEAST(Y, label, k, RFE.Filter.FEAST);

    case 2 % MRMR
        ind = nk_MRMR(Y, label, k, RFE.Filter.MRMR);
        if VERBOSE, fprintf(' MRMR'); end

    case 3 % Pearson
        if RFE.Filter.Pearson == 1, meas = 'pearson'; else, meas = 'spearman'; end
        weights = abs(nk_CorrMat(Y, label, meas));
        [~, ind] = sort(weights,'descend');
        if VERBOSE, fprintf(' %s', meas); end

    case 4 % Simba
        switch RFE.Filter.simba.utilfunc
            case 1 % linear
                if VERBOSE, fprintf(' Simba linear'); end
            case 2 % sigmoid
                if VERBOSE, fprintf(' Simba sigmoid (beta=%g)',RFE.Filter.simba.extra_param.beta); end
            case 3 % sigmoid with beta auto-aggust
                sugBeta = suggestBeta(Y,label);
                if VERBOSE, fprintf(' Simba sigmoid (beta=%g)',sugBeta); end
                RFE.Filter.simba.extra_param.beta = sugBeta;
        end
        weights = nk_SimbaMain(Y,label,RFE.Filter.simba.extra_param);
        [~, ind] = sort(weights,'descend');

    case 5 % Gflip
        switch RFE.Filter.gflip.utilfunc
            case 1
                if VERBOSE, fprintf(' G-flip zero-one'); end
            case 2
                if VERBOSE, fprintf(' G-flip linear'); end
            case 3
                if VERBOSE, fprintf(' G-flip sigmoid (beta=%g)', RFE.Filter.gflip.extra_param.beta); end
        end
        [~, weights] = gflip(Y, label, RFE.Filter.gflip.extra_param);
        [~, ind] = sort(weights,'descend');

    case 6 % AMS
        if VERBOSE, fprintf(' AMS'); end
        [weights, C] = ams(Y,label,RFE.Filter.ams.method, RFE.Filter.ams.sphere, ...
            RFE.Filter.ams.extra_param);
        [ ~, ind] = sort(weights,'descend');

    case 7 % IMRelief
        if iscell(RFE.Filter.imrelief)
            if numel(RFE.Filter.imrelief)<curclass
                imrelief = RFE.Filter.imrelief{end};
            else
                imrelief = RFE.Filter.imrelief{curclass};
            end
        else
            imrelief = RFE.Filter.imrelief;
        end
        nsigma = numel(imrelief.sigma);
        nlambda = numel(imrelief.lambda);
        if nsigma > 1 || nlambda > 1, weights = zeros(size(Y,2),nsigma,nlambda); end
        for i=1:nsigma
            for j=1:nlambda
                if VERBOSE, fprintf(' IMRelief (sigma=%g, lambda=%g)',  imrelief.sigma(i), ...
                                                            imrelief.lambda(j)); end
                weights(:,i,j) = IMRelief_Sigmoid_FastImple(Y', label,  ...
                                                    imrelief.distance, ...
                                                    imrelief.sigma(i), ...
                                                    imrelief.lambda(j), ...
                                                    imrelief.maxiter, ...
                                                    imrelief.plotfigure);
            end
        end
        % Perform crossvalidation
        if nsigma > 1 || nlambda > 1

            Mx = zeros(nsigma,nlambda);
            for i=1:nsigma
                for j=1:nlambda
                    pred = IMRelief_Sigmoid_FastImple_Predict2(Y', ...
                                                        label, ...
                                                        Ynew', ...
                                                        weights(:,i,j), ...
                                                        imrelief.distance, ...
                                                        imrelief.sigma(i));

                    Mx(i,j) = feval('BAC2', labelnew, pred);
                end
            end

            mxcv = max(Mx(:));
            [xpos, ypos] = find(Mx == mxcv);
            if VERBOSE, 
                fprintf('\nComputed CV grid:')
                for i=1:nsigma,
                    fprintf('\n')
                    fprintf('\t%1.1f',Mx(i,:)); 
                end
                fprintf('\nSelected sigma = %1.1f, lambda = %1.3f, Perf = %1.2f', ...
                    imrelief.sigma(xpos(1)), ... 
                    imrelief.lambda(ypos(1)), mxcv);
            end
            [~, ind] = sort(weights(:,xpos(1),ypos(1)),'descend');
            if exist('aMx','var')
                if isempty(aMx)
                    aMx = Mx;
                else
                    aMx = aMx + Mx;
                end    
            end
        else
            [~, ind] = sort(weights,'descend');
        end

    case 8

    case 9 % Increasing subspace cardinality

        ind = 1:k;
        if VERBOSE, fprintf(' increasing subspace cardinality (no filter)'); end

    case 10 % RGS (regression)
        if VERBOSE, fprintf(' RGS'); end
        %label = nk_PerfScaleObj(label);
        alpha = RGS(Y, label, RFE.Filter.RGS.extra_param);
        [~, ind] = sort(alpha,'descend');

    case 11 % Random subspace sampling
        if VERBOSE, fprintf(' Random subspace sampling'); end
        S = nk_RSS(size(Y,2), RFE.Filter.RSS);

    case 12
        if VERBOSE, fprintf(' RELIEF-based feature ranking'); end
        ind = relieff(Y,label,RFE.Filter.Relief.k);
        
    case 13
        weights = nk_FScoreFeatRank(Y, label);
        [~, ind] = sort(weights,'descend');
        if VERBOSE, fprintf(' FScore'); end
        
    case 14
        weights = nk_BhattacharyaFeatRank(Y, label);
         [~, ind] = sort(weights,'descend');
        if VERBOSE, fprintf(' Bhattacharya'); end

    
end