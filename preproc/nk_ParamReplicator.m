function [data_ind, train_ind, nP, nA] = nk_ParamReplicator(P, PXopt, PREPROC, TrainedParam)
% ==========================================================================
% [data_ind, train_ind] = nk_ParamReplicator(P, PXopt, PREPROC, TrainedParam)
% ==========================================================================
% nk_PerfParamMixerObj is a core NM function that allows to replicate
% TrainedParam structures according to overlap of the overall optimized 
% parameter array PXopt and the subset of PXopt consisting of the minimal 
% parameter combinations that were finally trained in the preprocessing
% pipeline
% ==========================================================================
% (c) Nikolaos Koutsouleris, 5/2016

% nP = Number of parameter combinations in PXopt, which may contain duplicates.
% Duplicates are due to unique ML parameter combinations associated with a
% single training data shelf.
nP = size(PXopt,1);

% nA = Number of preprocessing steps in pipeline
nA = numel(TrainedParam);

% nU = Number of unique parameter combinations that had to be trained in the
% preprocessing pipeline
nU = 1; emptfl = true;

% check whether parameters have to be replicated at all
if ~isempty(P); 
    Popt = P.opt;
    nU = size(Popt,1); 
    emptfl = false; 
end

% OK, we have to replicate
if ~emptfl
    
    indmat = zeros(size(Popt)); nPz = size(Popt,2);
    % Build index array to parameters
    for zz=nPz:-1:1
        uL = unique(Popt(:,zz));
        nuL = numel(uL);
        for zzz=1:nuL
            indmat(Popt(:,zz) == uL(zzz),zz) = zzz;
        end
    end
    
    tindmat = indmat;
    
    for zz=1:nPz-1
        uL1 = unique(indmat(:,zz));nuL1 = numel(uL1); cnt = 0;
        if zz>1
            uLx = unique(indmat(:,zz-1)); nuLx = numel(uLx); vec = [];
            for zzz=1:nuL1
                vec = [vec; repmat(uL1(zzz),nuLx,1)];
            end
        else
            vec = uL1;
        end
        summer = uL1(end); 
        for zzz=1:numel(vec):nU
            tindmat(zzz: zzz+numel(vec)-1,zz) = vec + summer*cnt;
            cnt = cnt+1;
        end
    end
    
    indmat = ones(nU, nA); fnd = false; cnt = nPz; zz=1;
    prevec = ones(nU,1);
    while zz <= nA
        if cnt > 0
            if strcmp(PREPROC.ACTPARAM{zz}.cmd, P.cmd{cnt}) || strcmp(P.cmd{cnt},'spatialfilter'), 
                prevec = tindmat(:,cnt); cnt = cnt-1;
            end
        end
        indmat(:,zz) = prevec;    
        zz=zz+1;
    end
end
data_ind = zeros(nP,1);
train_ind = zeros(nP,nA);

% Find pointers of current parameter mixer rows to PXopt
% parameter rows
if nP == 1 || emptfl
    data_ind(1:nP) = 1;
    train_ind(1:nP,1:nA) = 1;
else
    if size(PXopt,2) > size(Popt,2), 
        PXopt = rem_unique_cols(PXopt);
    end
    if size(PXopt,2) < size(Popt,2), 
        Popt = rem_unique_cols(Popt);
    end
    for z = 1:nU
        try
            ind_z = ismember(PXopt, Popt(z,:),'rows'); 
        catch
            fprintf('error')
        end
        if any(ind_z) 
            data_ind(ind_z) = z;
            train_ind(ind_z,:) = repmat(indmat(z,:),sum(ind_z),1);
        end
    end                    
end

end

function Pm = rem_unique_cols(Pm)

   remvec=[];
   for z = 1:size(Pm,2)
       Pm_c = Pm(:,z);
       if numel(unique(Pm_c))==1
           remvec=[remvec z];
       end
   end
   Pm(:,remvec)=[];
   
end
