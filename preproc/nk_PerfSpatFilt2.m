function sY = nk_PerfSpatFilt2( Y, CURACT, Param )

if isfield(CURACT,'SPATIAL') && CURACT.SPATIAL.cubetype>1
    % If data fusion has happened before the smoothing
    % extract the modality, smooth it and fuse it again 
    S = CURACT.SPATIAL;
    dimvecx = Param.dimvecx;
    if iscell(Param.brainmask), nM = numel(Param.brainmask); else, nM=1; end
    nP = size(CURACT.SPATIAL.PX.opt,1);
    if ~nP, nP=1;end
    fY = cell(nP,nM);
    sY = cell(nP,1);
    Vm = []; Vmvol = [];
    for zu = 1:nM
        if nM>1
            brainmask       = nk_RemCommaFromStr(Param.brainmask{zu});
            datatype        = Param.datatype(zu);
            if isfield(Param,'Vm')
                Vm          = Param.Vm{zu};
                Vmvol       = Param.Vmvol{zu};
            end
            badcoords       = Param.badcoords{zu};
            LABEL           = Param.threshval{zu};
            LABELOPERATOR   = Param.threshop{zu};
        else
            brainmask       = nk_RemCommaFromStr(Param.brainmask{1});
            datatype        = Param.datatype;
            if isfield(Param,'Vm')
                Vm          = Param.Vm;
                Vmvol       = Param.Vmvol;
            end
            badcoords       = Param.badcoords{1};
            LABEL           = Param.threshval;
            LABELOPERATOR   = Param.threshop;
        end
        dimst               = dimvecx( zu ) + 1; 
        dimend              = dimvecx( zu + 1 );
        tY                  = Y(:,dimst:dimend);     
        if datatype ~= 1
            for zv = 1:nP
                fY{zv, zu} = tY; 
            end
        else
            tmpflg = false;
            if ~exist(brainmask,'file')
                if ~isempty(Vm) && ~isempty(Vmvol)
                    tmpflg = true;
                    Vm = WriteTempVol(Vm, Vmvol, tmpflg);
                    brainmask = Vm.fname;
                else
                    error('Space-defining image %s cannot be found! Make sure your paths are up-to-date!', brainmask);
                end
            end
            S.brainmask                  = brainmask;
            S.badcoords                  = badcoords;
            S.Vm                         = spm_vol(S.brainmask);
            [S.dims, S.indvol, ~, S.vox] = nk_ReadMaskIndVol(S.Vm, [], LABEL, LABELOPERATOR);
            if isfield(Param,'indNonRem') && ~isempty(Param.indNonRem) && sum(~Param.indNonRem) > 0
                indNonRem                = Param.indNonRem;
                ttY                      = zeros(size(tY,1),size(indNonRem)); 
                ttY(:,indNonRem)         = tY; 
            else
                ttY                      = tY;
            end
            nanMask                      = isnan(ttY);
            indnan                       = any(nanMask,2);
            tY                           = nk_SpatFilt(ttY, S);
            clear ttY
            % Transfer smoothed data into container and handle NaNs
            % properly
            for zv = 1:nP, 
                fY{zv,zu} = tY{zv};
                if any(indnan), 
                    fY{zv,zu}(nanMask) = nan; 
                end
            end
            if tmpflg, DeleteTempVol(Vm,tmpflg); end 
        end
    end
    for zv = 1:nP, sY{zv} = cell2mat(fY(zv,:));end
else
    sY = Y;
end