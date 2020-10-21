function PA = nk_PermX(inp, strout, id, GridAct, batchflag)

global SVM RFE SAV DR GRD MULTI RAND SCALE DISCRET SYMBOL COVAR TEST MODEFL

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%% INITIALIZATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
nclass          = inp.nclass;       % # of binary comparisons
label          = inp.label;
lx              = inp.l;            % # of subjects
cv              = inp.cv;           % cross-validation structure
dimension       = inp.dim_index;    % current dimension to analyze
preprocmat      = inp.preprocmat;   % paths to preprocessed data files
gdmat           = inp.gdmat;        % paths to precomputed GD structures
nperms          = inp.nperms;       % # of permutations
ovrwrt          = inp.ovrwrt;       % overwrite existing data files
analysis        = inp.analysis;     % GDanalysis structure to be used
% Check whether RVM or SVM will be used?
if strcmp(SVM.prog,'MikSVM'), algostr = 'RVM'; else algostr = 'SVM'; end

% Setup CV2 container variables:
[ix, jx] = size(cv.TrainInd);

if ~exist('GridAct','var') || isempty(GridAct), ...
        GridAct=nk_CVGridSelector(ix,jx); end

if ~exist('batchflag','var') || isempty(batchflag), batchflag = false; end

%%%%%%%%%%%%%%%%% SETUP (P)ermutation (A)nalysis STRUCTURE %%%%%%%%%%%%%%%%

if ~batchflag
    PA.params.SVM       = SVM;
    PA.params.RFE       = RFE;
    PA.params.SAV       = SAV;
    PA.params.GRD       = GRD;
    PA.params.COVAR     = COVAR;
    PA.params.SCALE     = SCALE;
    PA.params.DISCRET   = DISCRET;
    PA.params.SYMBOL    = SYMBOL;
    PA.params.DR        = DR;
    PA.params.TEST      = TEST;
    PA.params.MULTI     = MULTI;
    PA.params.RAND      = RAND;
    PA.perms            = nperms;
    PA.dimension        = inp.dimension;
    PA.dimension_index  = inp.dim_index;
    PA.GridAct          = GridAct;
    PA.GDpaths          = cell(ix,jx);
    PA.SvmBinComp       = nclass;
end

ll = 1; ol = 0;
signTs = nan(ix, jx, nclass); Ts = nan(ix, jx, nclass); 
Nx = nan(ix, jx, nclass); Perms = nan(ix, jx, nclass);

if strcmp(MODEFL,'classification')
    o_pred = cell(lx,nclass);
    m_pred = cell(lx,nclass,nperms);
    binmode = 1;
    if MULTI.flag, 
        signmTs = zeros(ix,jx);
        omulti_pred = cell(lx,1);
        pmulti_pred = cell(lx,nperms);
    end;
elseif strcmp(MODEFL,'regression')
    o_pred = cell(lx,1);
    m_pred = cell(lx,nperms);
    binmode = 0;
end

switch SVM.GridParam
    case {9,11,12}
        act = '';
    otherwise
        act = 'inv';
end
for f=1:ix % Loop through CV2 permutations

    for d=1:jx % Loop through CV2 folds
        
        if ~GridAct(f,d), ll = ll + 1; continue, end;
        TsInd = cv.TestInd{f,d};
        operm = f; ofold = d;
        
        %%%%%%%%%%%%%%%%%%%%%%%%% USE PRECOMPUTED ? %%%%%%%%%%%%%%%%%%%%%%%
        cvstr = ['_oCV' num2str(f) '.' num2str(d) '_'];
        oPERMpath = fullfile(pwd,[SAV.matname strout '_PERMdatamat' cvstr 'ID' id '.mat']);
        if exist(oPERMpath,'file') && ~ovrwrt && ~batchflag
            [opth, onam] = fileparts(oPERMpath);
            fprintf('\nPERMdatamat found for CV2 [%g,%g]:',f,d)
            fprintf('\nLoading: %s',onam)
            load(oPERMpath)
            signTs(f,d,:) = signTsh; Nx(f,d,:) = Nxh; Perms(f,d,:) = Permsh;
            %ECDf = ECDf+ECDfh;
            if MULTI.flag, signMult(f,d) = Mult; end
            continue
        elseif exist(oPERMpath,'file') && ~ovrwrt && batchflag
            [opth, onam] = fileparts(oPERMpath);
            fprintf('\nPERMdatamat found for CV2 [%g,%g]:\n%s',f,d,onam)
            fprintf('\nBatch mode detected. Continue.')
            continue
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%% LOAD PREPROCDATA %%%%%%%%%%%%%%%%%%%%%%%%%
        pppth = deblank(preprocmat{1,f,d});
        [ppth,pnam] = fileparts(pppth);
        if isempty(pppth) || ~exist(pppth,'file')
            warning(['No valid PreprocData-MAT detected for CV2 partition ' ...
                '[' num2str(f) ', ' num2str(d) ']!']);
            continue
        else
            fprintf('\n\nLoading preprocessed data:');
            fprintf('\n%s',pnam);
            load(pppth,'mapY')
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%% LOAD CVDATAMAT  %%%%%%%%%%%%%%%%%%%%%%%%
        gdpth = deblank(gdmat{1,f,d});
        if isempty(gdpth) || ~exist(gdpth,'file')
            warning(['No valid CVdatamat-MAT detected for CV2 partition ' ...
                '[' num2str(f) ', ' num2str(d) ']!']);
            continue
        else
            [gpth,gnam] = fileparts(gdpth);
            fprintf('\n\nLoading corresponding CVdatamat:');
            fprintf('\n%s',gnam);
            load(gdpth)
        end
        
        ol = ol+1;
        [iy, jy] = size(cv.cvin{f,d}.TrainInd);
        
        %%%%%% Build permutation structure for current CV2 partition %%%%%%
        oGD     = cell(nclass,1); pGD = cell(nclass,nperms);
        if RFE.ClassRetrain
            llx     = size(cv.cvin{f,d}.TrainInd{1,1},1) + size(cv.cvin{f,d}.TestInd{1,1},1) ;
        else
            llx     = size(cv.cvin{f,d}.TrainInd{1,1},1);
        end
        pTrInd  = zeros(llx,nperms);
        
        fprintf('\n-------------------------------------------------------------------------\n')
        % Permute m times the order of cases in current CV2 training partition
        fprintf('\nBuild permutation index structure (nperms=%g)\n',nperms); 
        for m=1:nperms
            pTrInd(:,m) = randperm(llx); 
        end
        
        for h=1:nclass % Loop through binary comparisons
            
            if nclass > 1, fprintf('\n\n*** %s #%g ***',algostr, h); end

            switch binmode
                case 0 % MULTI-GROUP or regression analysis
                    TsI = TsInd;
                case 1 % BINARY-GROUP analysis
                    TsI = TsInd(cv.classnew{f,d}{h}.ind);
                    if min(mapY{h}.TsL)~=-1 || max(mapY{h}.TsL) ~= 1
                        mapY{h}.TsL = cv.classnew{f,d}{h}.label;
                    end      
            end
            
            if MULTI.flag
                % Get optimum binary classifier params from GDanalysis
                if inp.nocexp
                    C = NaN; 
                else
                    m = size(analysis.multi_bestcpos,1);
                    if m > 1
                        C = analysis.C{h}(ll,analysis.multi_bestcpos(ll)); 
                    else
                        C = analysis.C{h}(analysis.multi_bestcpos(ll)); 
                    end
                end
                if inp.nogexp
                    G = NaN; 
                else
                    m = size(analysis.multi_bestgpos,1);
                    if m > 1
                        G = analysis.Gamma{h}(ll,analysis.multi_bestgpos(ll)); 
                    else
                        G = analysis.Gamma{h}(analysis.multi_bestgpos(ll)); 
                    end
                end

            else
                % Get optimum binary classifier parameters from GDanalysis
                if inp.nocexp,
                    C = NaN; 
                else
                    m = size(analysis.bestc{h},1);
                    if m > 1
                        C = analysis.bestc{h}(ll);  
                    else
                        C = analysis.bestc{h};  
                    end
                end
                % Get optimum binary classifier parameters from GDanalysis
                if inp.nocexp,
                    G = NaN; 
                else
                    m = size(analysis.bestg{h},1);
                    if m > 1
                        G = analysis.bestg{h}(ll);  
                    else
                        G = analysis.bestg{h};  
                    end
                end   
            end
            
            if ~strcmp(SVM.prog,'CUDSVM') && ~strcmp(SVM.prog,'MikRVM')
                kstr = num2str(C,'%1.10f');
                if isnan(G), gstr='0'; else gstr = num2str(G,'%1.10f'); end
            else
                kstr = C; if isnan(G), gstr = 0; else gstr = G; end
            end
            
            MD = cell(iy,jy);
            
            if MULTI.flag
                oGD{h}.DT{1,1} = cell(iy,jy); oGD{h}.RT{1,1} = cell(iy,jy);
                oGD{h}.DV{1,1} = cell(iy,jy); oGD{h}.RV{1,1} = cell(iy,jy);
                
            end
            
            
            %%%%%%%%%%%%%%%%%%%%%% MOVE TO CV1-LEVEL %%%%%%%%%%%%%%%%%%%%%%
            for k=1:iy

                for l=1:jy
                    
                    % Get Feature Subspace Indices
                    ul = size(GD{h}.bestfeats{k,l},2);
                    
                    if ul > 1
                        fprintf('\nCV2 [%g,%g], CV1 [%g,%g], Compute original models (ensemble model N: %g) ...', f,d,k,l, ul)
                    else
                        fprintf('\nCV2 [%g,%g], CV1 [%g,%g], Compute original models (single model) ...', f,d,k,l)
                    end
                    
                    MD{k,l} = cell(ul,1); 
                    oGD{h}.DT{k,l} = zeros(llx,ul);
                    oGD{h}.DV{k,l} = zeros(size(mapY{h}.CV{k,l},1),ul);
                    
                    %%%%%%%%%%%%% RECONSTRUCT ORIGINAL MODEL %%%%%%%%%%%%%%
                    % Reconstruct original classifier (ensemble)
                    for u=1:ul
                        
                        if ~islogical(GD{h}.bestfeats{k,l}(:,u))
                            F = GD{h}.bestfeats{k,l}(:,u) ~= 0;
                        else
                            F = GD{h}.bestfeats{k,l}(:,u);
                        end
                        
                        XTr     = mapY{h}.Tr{k,l}(:,F);
                        XCV     = mapY{h}.CV{k,l}(:,F);
                        TrL     = mapY{h}.TrL{k,l};
                        CVL     = mapY{h}.CVL{k,l};
                        
                        if RFE.ClassRetrain, XTr = [XTr;  XCV]; TrL = [TrL; CVL]; end
                        
                        % Model computation:
                        [perf, MD{k,l}{u}] = nk_GetParam(XTr, TrL, kstr, gstr);

                        if MULTI.flag
                            oGD{h}.DT{1,1}{k,l} = perf.dec_values; oGD{h}.RT{1,1}{k,l} = perf.target;
                            [dum, oGD{h}.DV{1,1}{k,l}, oGD{h}.RV{1,1}{k,l}] = ...
                                nk_GetTestParam(XTr, [], XCV, CVL, MD{k,l});
                        end
                        
                    end
                end
                
            end
            %%%%%%%%%%%%% RECOMPUTE ORIGINAL CV2-PERFORMANCE %%%%%%%%%%%%%%
            
            % Check if weight structure exists
            if ~isfield(GD{h},'bestweights')
                Weights = cell(iy, jy);
                for k=1:iy
                    for l=1:jy
                        ul = size(GD{h}.bestfeats{k,l},2);
                        Weights{k,l} = ones(ul,1);
                    end
                end 
            else
                Weights = GD{h}.bestweights;
            end
            
            fprintf('\n\nCV2 [%g,%g], Evaluate CV2 performance of original models ... ', f, d)
            
            [TsP, ds, rs, hrs, hds, hrx, hdx, hTsP_rs, hTsP_ds] = ...
            nk_GetTestParam(mapY{h}.Tr, mapY{h}.CV, mapY{h}.Ts, mapY{h}.TsL, MD, GD{h}.bestfeats, Weights);
            
            if MULTI.flag, oGD{h}.DS{1,1} = ds; oGD{h}.RS{1,1} = rs; end
            % Compute mean CV2 performance across CV1 grid
            mTsP = mean(nk_cellcat(TsP));
        
            %%%%%%%%%%%%%%%%% COMPUTED PERMUTED MODELS %%%%%%%%%%%%%%%%%%%%           
            fprintf('\n\nPermutation analysis in CV2 [%g,%g]: %g perms', f, d,nperms)
            
            % Allocate memory:
            phTsP_rs    = zeros(nperms,1); 
            phTsP_ds    = zeros(nperms,1);
            mpTsP       = zeros(nperms,1);
            
            for m=1:nperms
                
                pMD = cell(iy,jy);
                
                for k=1:iy
                
                    for l=1:jy
%                           
                        % Get Feature Subspace Indices
                        ul = size(GD{h}.bestfeats{k,l},2);

                        if ul > 1
                            fprintf(['\n%g / %g perms: CV2 [%g,%g], CV1 [%g,%g],' ...
                                ' Compute permuted models (ensemble N: %g) ...'],m, ...
                                nperms, f,d,k,l,ul)
                        else
                            fprintf(['\n%g / %g perms: CV2 [%g,%g], CV1 [%g,%g],' ...
                                ' Compute permuted models (single classifier) ...'], m, ...
                                nperms, f,d,k,l)
                        end
                        % Reconstruct original classifier (ensemble)
                        for u=1:ul
                            pMD{k,l} = cell(ul,1);
                            pGD{h}.DT{k,l} = zeros(llx,ul);
                            pGD{h}.DV{k,l} = zeros(size(mapY{h}.CV{k,l},1),ul);

                            % Permuted model computation:
                            [perf, pMD{k,l}{u}] = nk_GetParam(XTr(pTrInd(:,m),:), TrL, kstr, gstr);

                            if MULTI.flag
                                pGD{h,m}.DT{1,1}{k,l} = perf.dec_values; pGD{h,m}.RT{1,1}{k,l} = perf.target;
                                [dum, pGD{h,m}.DV{1,1}{k,l}, pGD{h,m}.RV{1,1}{k,l}] = ...
                                    nk_GetTestParam(XTr(permind,:), [], XCV, CVL, pMD{k,l});
                            end
                        end
                    end
                end
                fprintf('\n\nCV2 [%g,%g], Evaluate CV2 performance of permuted models ... ', f, d)
                
                [pTsP, pds, prs, phrs, phds, phrx, phdx, phTsP_rs(m), phTsP_ds(m)] = ...
                nk_GetTestParam(mapY{h}.Tr, mapY{h}.CV, mapY{h}.Ts, mapY{h}.TsL, pMD(:,:,m), GD{h}.bestfeats, Weights);
               
                if MULTI.flag, pGD{h,m}.DS{1,1} = pds; pGD{h,m}.RS{1,1} = prs; end
                % Eventually, compute mean permuted CV2 performance across
                % CV1 grid:
                if ~RFE.CV2Class.EnsembleStrategy.Metric, mpTsP(m) = nk_cellcat(pTsP); end
            
            end
            
            fprintf(' done.\n')
            %%%%%%%%%%% COMPUTE SIGNIFICANCE (OF BIN. PROBLEM) %%%%%%%%%%%% 
            % Choose performance type:
            switch RFE.CV2Class.EnsembleStrategy.Metric
                case 0 % Mean CV2 performance across CV1 grid
                    [signTs(f,d,h), Nx(f,d,h)] = ...
                        compute_probability(mpTsP, mTsP, act);
                    Ts(f,d,h)     = mTsP;
                case 1 % Ensemble performance (target values)
                    [signTs(f,d,h), Nx(f,d,h)] = ...
                        compute_probability(phTsP_rs, hTsP_rs, act);
                    Ts(f,d,h)     = hTsP_rs;
                case 2 % Ensemble performance (decision values)
                    [signTs(f,d,h), Nx(f,d,h) ] = ...
                        compute_probability(phTsP_ds, hTsP_ds,act);
                    Ts(f,d,h)     = hTsP_ds;
            end
            %ECDf(:,h) = ECDf(:,h) + ecdfh;
            Perms(f,d,h) = nperms;
        end
        
        %%%%%%%%%%% COMPUTE SIGNIFICANCE OF MULTI-CLASS PROBLEM %%%%%%%%%%%
        if MULTI.flag
            
            % Original multi-class model:
            M = nk_MultiClassPermFold(oGD, f, d, false, nclass);
            omTs = M.mTsPerf;
            omulti_pred(TsInd) = cellmat_mergecols(omulti_pred(TsInd), ...
                                num2cell(M.TsPred{1,1},2));
            pmTs = zeros(nperms,1);
            
            for m=1:nperms
                
                M = nk_MultiClassPermFold(pGD(:,m), f, d, false, nclass);
                pmTs(m) = M.mTsPerf;
                pmulti_pred(TsInd,m) = cellmat_mergecols(pmulti_pred(TsInd,m), ...
                                num2cell(M.TsPred{1,1},2));
            end
                    
            [signmTs(f,d), mNx(f,d) ] = ...
                        compute_probability(pmTs, omTs, act);
        end
        
        %%%%%%% SAVE PERMUTATION ANALYSIS FOR CURRENT CV2 PARITION %%%%%%%%
        %if batchflag
            fprintf('\nSaving %s', oPERMpath)
            Nxh = Nx(f,d,:); signTsh = signTs(f,d,:); Permsh = Perms(f,d,:);
            %ECDfh = ECDf./ol;
            if MULTI.flag
                TsM = signmTs(f,d);
                save(oPERMpath,'signTsh','Nxh','Permsh', ...
                    'TsM', 'kstr','gstr','operm','ofold');
            else
                save(oPERMpath,'signTsh','Nxh','Permsh','kstr','gstr','operm','ofold');
            end
        %end
        ll=ll+1;
    end
end

if ~batchflag
    for h=1:nclass
        PA.TS{h}            = Ts(:,:,h);
        PA.sign_TS{h}       = signTs(:,:,h);
        
        PA.mean_sign_TS{h}  = mean(PA.sign_TS{h}(~isnan(PA.sign_TS{h})));
        PA.std_sign_TS{h}   = std(PA.sign_TS{h}(~isnan(PA.sign_TS{h})));
        PA.omni_sign_TS{h}  = Nx(:,:,h); SP = Perms(:,:,h);
        PA.omni_sign_TS{h}  = sum(PA.omni_sign_TS{h}(~isnan(PA.omni_sign_TS{h}))) / sum(SP(~isnan(SP)));
        %PA.ECDf{h}          = ECDf(:,h)./ol;
    end
    O = zeros(lx,1); P = zeros(lx,nperms); Perr = zeros(nperms,1);
    if MULTI.flag
        for i=1:lx
            % Compute mean across CV2 perms of original multi-class model
            if isempty(omulti_pred{i})
                O(i) = NaN; 
            else
                O(i) = mean(omulti_pred{i});
            end
            % Compute mean across CV2 perms of permuted multi-class models
            for j=1:nperms
                if isempty(pmulti_pred{i,j})
                    P(i,j) = NaN;
                else
                    P(i,j) = mean(pmulti_pred{i,j});
                end
            end
            P = round(P); O = round(O);
            Oerr = (lx-sum(label~= O))*100;
            for j=1:nperms
                Perr(j) = (lx-sum(label~= P(:,j)))*100;
            end
            MTS = compute_probability(Perr, Oerr, stattool, act);
        end
    end
else
    PA = [];
end

return

% ________________________________________________________________________
function [Px, Nx] = compute_probability(perms, obs, act)

if strcmp(act,'inv')
    Nx = sum(perms >= obs);
else
    Nx = sum(perms <= obs);
end
Px = (Nx / length(perms));

return
