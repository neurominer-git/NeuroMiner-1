function nk_FeatSel_Cluster(dat, xs)
%
% function nk_FeatSel_Cluster(dat)
% ================================

% (c) 2009, Nikolaos Koutsouleris, last modified 02/10/2009

[ix,jx] = size(dat.cv.TrainInd);
brainmask = regexprep(dat.brainmask,',1','');
cv = dat.cv;
id = dat.id;

if strcmp(spm('ver'),'SPM5'),  
    spm5ver = 1; 
else
    spm5ver=0;
end

if ~exist(brainmask,'file')
	
	if spm5ver
        brainmask = spm_select(1,'image','select brainmask');
    else
        brainmask = spm_get(1,'IMAGE','select brainmask');
    end

end

Vm = spm_vol(brainmask);
[dims, indvol] = nk_ReadMaskIndVol(Vm);

[featmat,FeatInfo] = load_mats(ix,jx, spm5ver, id);

% Retrieve FEATSEL
load(featmat{1},'FEATSEL');

% Select image volumes to write
imgvec = zeros(length(FEAT{1,1}.dat),1);
for i = 1:length(FEAT{1,1}.dat)
	str1 = FEAT{1,1}.dat{i}.datstr;
	str2 = FEAT{1,1}.dat{i}.datdesc;
	fprintf('Image volume %g: %s (%s)',i, str1, str2);
end
imgind = nk_input('Select image volumes to write to disk',0,'e');
imgvec(imgind) = 1;

clc;
clear dat

% Build output data string
switch FEATSEL.symflag
    case 1
        symsuff = '_dsc';
    case 2
        symsuff = '_sym';
    otherwise
        symsuff = '';
end

switch FEATSEL.varflag
    case 1
        switch FEATSEL.ranktype 
            case 1, 
                strtype = ['pearson' symsuff];
            case 2
                strtype = ['FDR' symsuff];
            case 3
                strtype = ['FTest' symsuff];
            case 4
                strtype = 'MI';  
        end
    case 2
        switch FEATSEL.ranktype 
            case 1, 
                strtype = ['PLS' symsuff];
            case 2
                strtype = ['kPLS' symsuff];
            case 3
                strtype = ['PCA' symsuff];
        end
end

switch FEATSEL.cubetype
    case 1
        strtype = [strtype '_sp4x4'];
    case 2
        strtype = [strtype '_sp9x3'];
end

fprintf('\n')
fprintf('\n\nCluster Extraction')
fprintf('\n------------------')
fprintf('\n')
%fprintf('\nParams:')
FEATSEL
switch FEATSEL.varflag
    case 1
        if ~FEATSEL.permflag 
            permsuff = '';
        else
            permsuff = ['_p' num2str(FEATSEL.nperms)];
        end
        switch FEATSEL.cv,
            case 0
                suff = permsuff;
            case 1
                suff = ['_loo' permsuff];
            case 2
                suff = ['_cv' num2str(FEATSEL.kfold) permsuff];
            case 3
                suff = ['_bs' num2str(FEATSEL.nboot) permsuff ];
        end
    case 2
        suff = [];
        if FEATSEL.bootflag 
            suff = ['_bs' num2str(FEATSEL.nboot)];
        end
        if FEATSEL.permflag
            suff = [suff '_p' num2str(FEATSEL.nperms)];
        end
end

Vfeat = Vm;
[iy,jy] = size(cv.cvin{1,1}.TrainInd);
if nargin < 2
  xs = continue_interrupted;
  %[xs, ys, FEAT, FEATSEL] = continue_interrupted(spm5ver, id, FEATSEL);
  opermvec=[xs(1):ix];
  ofoldvec=[xs(2):jx];
  ipermvec=1:iy;
  ifoldvec=1:jy;
else
% this is for batch mode
  %ys=1;
  opermvec=xs(1);
  ofoldvec=xs(2);
  ipermvec=1:iy;
  ifoldvec=1:jy;
end

for iz=1:length(opermvec) % Loop through all outer CV samples [i=perm,j=fold]

	for jz=1:length(ofoldvec)
		
		i = opermvec(iz);
		j = ofoldvec(jz);
        
		for v=1:length(FeatInfo)
			if FeatInfo{v}.CVInd(1) == iz && FeatInfo{v}.CVInd(2) == jz
				fprintf('\n\nLoading %s.',FeatInfo{v}.matname);
				load(deblank(featmat{v})); % load FEAT information for current outer perm / outer fold
				break
			end
		end
        
		for kz=1:length(ipermvec) % Loop through all inner CV samples [i,j],[k=inner perm,l=inner fold]
            
			for lz=1:length(ifoldvec)

				k = ipermvec(kz);
				l = ifoldvec(lz);
       
                fprintf('\n\nWriting images & extracting clusters of CV2 (perm,fold=[%g,%g]), CV1 (perm, fold=[%g,%g]).',i,j,k,l)
                cvstr = ['oCV' num2str(i) '.' num2str(j) '_iCV-' num2str(k) '.' num2str(l)];

                Vfeat.dt=[4,0]; % save in INT16
                
                xu = size(FEAT{k,l}.dat{1}.vol,2);
                
                if FEATSEL.writeimg

                    for u=1:xu
                        if xu>1, 
                            pref = [FEATSEL.prefix '_cl' num2str(u) '_' strtype];
                        else
                            pref = [FEATSEL.prefix '_' strtype];
                        end
                        cntz=1;
                        for z=1:length(FEAT{k,l}.dat)
			   if imgvec(z) 
				if isfield(FEAT{k,l}.dat{z},'tag'),
					if ~strcmp(FEAT{k,l}.dat{z}.tag,'clust')
						dat{cntz} = FEAT{k,l}.dat{z};
						cntz=cntz+1;
					else
						continue
					end
				else
					dat{cntz} = FEAT{k,l}.dat{z};
					cntz=cntz+1;
				end;
				
				if iscell(FEAT{k,l}.dat{z}.vol)
					vol=FEAT{k,l}.dat{z}.vol{u};
				else
					vol=FEAT{k,l}.dat{z}.vol(:,u);
				end
				
				FEAT{k,l}.dat{z}.fname{u,1} = writeimg(Vfeat, ...
					vol, ...
					FEAT{k,l}.dat{z}.datstr, ...
					FEAT{k,l}.dat{z}.datdesc, pref, cvstr, suff, dims, indvol, 1);
			   else
				FEAT{k,l}.dat{z}.fname{u,1} = fullfile(pwd, [pref '_' FEAT{k,l}.dat{z}.datstr '_' cvstr suff '.img']);
		           end
                        end
                        FEAT{k,l}.dat = dat;
                    end
                end
                
                % Use the written voxel-by-voxel discriminative volume to
                % generate a cluster-labeled map. It is highly advisable to
                % apply a spatial consistency criterion to the
                % voxel-by-voxel discriminative maps.
                datl = length(FEAT{k,l}.dat)+1;
                if FEATSEL.clustflag
                    for z=1:length(FEAT{k,l}.dat) % Get the sig volume
                        if strcmp(FEAT{k,l}.dat{z}.tag,'sig')
                            
                            % This routine returns 3D-volumes containing
                            % the cluster discriminative metric for each
                            % binary comparison
                            
                            clust = nk_Cluster(FEAT{k,l}.dat{z}.fname, cv.class{i,j}{u}.groupdesc);
                            
                            % Reshape clust into a vector;
                            % clustvec = zeros(ix,1);
                            clustvec = cell(xu,1);
                            for u=1:xu
                                
                                clustvec{u} = clust{u}(indvol);
                                if xu>1, 
                                    pref = [FEATSEL.prefix '_cl' num2str(u) '_' strtype];
                                else
                                    pref = [FEATSEL.prefix '_' strtype];
                                end

                                FEAT{k,l}.dat{datl}.vol{u} = clustvec{u};
                                FEAT{k,l}.dat{datl}.tag = 'clust';
                                FEAT{k,l}.dat{datl}.datstr = 'clust';
                                FEAT{k,l}.dat{datl}.datdesc = 'Clustered discriminative criterion';
                                FEAT{k,l}.dat{datl}.fname{u} = writeimg(Vfeat, ...
                                    FEAT{k,l}.dat{datl}.vol{u}, ...
                                    FEAT{k,l}.dat{datl}.datstr, ...
                                    FEAT{k,l}.dat{datl}.datdesc, pref, cvstr, suff, dims, indvol);
                            end
                            break;
                        end
                       
                    end
                end
            end
        end
        
        % Inner CV loop is finished here, (optionally) save results to disk
        cvstr = ['_oCV' num2str(i) '.' num2str(j)];outfold = i; outperm=j;
        matname = [FEATSEL.matname '_' strtype cvstr suff '_ID' id '.mat'];
        matpath = fullfile(pwd,matname);
        fprintf('\n\nSaving feature selection results of outer CV partition [fold = %g, perm = %g].',outfold,outperm);
        fprintf('\n==> %s',matname);
        save(matpath, 'id', 'FEATSEL','FEAT','outfold','outperm');
        clear FEAT; %ys=[1,1]; % Reset FEAT and counters
	end
end

% Save cross-validation settings in case it is necessary to reconstruct the
% analysis prerequisites
matname = [FEATSEL.matname '_' strtype '_CVpartData_ID' id '.mat'];
matpath = fullfile(pwd, matname);
fprintf('\n\n Saving cross-validation definitions to disk.')
fprintf('\n==> %s',matname);
save(matpath,'id','cv');
fprintf('\n\nDone.')

return

% _________________________________________________________________________
function fpath = writeimg(Vfeat, vol, datstr, datdesc, pref, cvstr, suff, dims, indvol, scaleflag)

if nargin < 10
    scaleflag=0;
end
fpath = fullfile(pwd, [pref '_' datstr '_' cvstr suff '.img']);
[fp,fn] = fileparts(fpath); 
Vfeat.fname = fpath;
Vfeat.descrip = datdesc;
if scaleflag
    Vfeat.pinfo(1) = ... & adjust scaling
        (max(vol(:,1)) - min(vol(:,1)))/32768;
end
Vdat = zeros(dims);
Vdat(indvol) = vol(:,1);
fprintf('\nWriting %s.', fn)
spm_write_vol(Vfeat,Vdat);

return


% ____________________________________________________
function xs = continue_interrupted
%function [xs, ys, FEAT, FEATSEL]= continue_interrupted(spm5ver, datid, FEATSEL)

xs=[1,1]; ys=[1,1]; FEAT=[];
 
xsflag = nk_input('Start from a specific CV2 (outer) partition?',0,'yes|no',[1,0],0);

if xsflag
    
    xs = nk_input('Indicate CV2 start position [outer perm, outer fold]',0,'e',[],2);
%      ysflag = nk_input('Start from a specific CV1 (inner) partition?',0,'yes|no',[1,0],0);
%      
%      if ysflag
%          ys = nk_input('Indicate CV1 start position [inner perm, inner fold]',0,'e',[],2);
%          loadflag = nk_input('Continue previous analysis?',0,'yes|no',[1,0],0);
%  
%          if loadflag
%              if spm5ver
%                  featmat=spm_select(1, ...
%                      ['oCV. ' num2str(xs(1)) '.' num2str(xs(2)) '.*ID' datid '\.mat'], ...
%                      'Select feature MATs');
%              else
%                  featmat=spm_get(1, ...
%                      ['*oCV' num2str(xs(1)) '.' num2str(xs(2)) '*ID' datid '.mat'], ...
%                      'Select feature MAT');
%              end
%              fprintf('\n\tLoading %s.',featmat);
%              load(featmat); 
%          end
%      end
end

return

function [featmat, FeatInfo] = load_mats(ix,jx, spm5ver, datid)

matnum = ix*jx;
if spm5ver
    featmat=spm_select(matnum,'oCV.*ID.*\.mat',['Select ' num2str(matnum) ' feature MATs']);
else
    featmat=spm_get(matnum,'*oCV*ID*.mat',['Select ' num2str(matnum) ' feature MATs']);
end

featmat = cellstr(featmat);
FeatInfo = cell(size(featmat,1),1);

for v=1:length(featmat)
    [pth,nam,ext] = fileparts(deblank(featmat{v}));
    FeatInfo{v}.id = regexp([nam ext],'ID.*\.','match');
    stoCV = regexp(nam,'oCV\d+\.\d+','match');

    if isempty(FeatInfo{v}.id)
        fprintf('\n\tNo ID found in filename %s', nam);
        fprintf('\nOpening mat file to find ID.');
        load(deblank(featmat{v}));
        try 
            FeatInfo{v}.id = id;
        catch
            error('\nNo ID found. Abort!')
        end
    end

    if isempty(stoCV)
        fprintf('\n\tNaming convention of mat file is unexpected.');
        fprintf('\n\tOpening mat file to find CV indices ...');
        load(deblank(featmat{v}));
        try 
            FeatInfo{v}.CVind(1) = outperm;
            FeatInfo{v}.CVind(2) = outfold;
        catch
            error('\n\tNo CV indices found. Abort!');
        end
    end

    FeatInfo{v}.id = char(FeatInfo{v}.id);
    FeatInfo{v}.id = FeatInfo{v}.id(3:end-1); % Extract ID number

    stoCV = char(stoCV);
    stoCV = stoCV(4:end);
    o = regexp(stoCV,'.\d');
    FeatInfo{v}.CVInd(1) = str2double(stoCV(1:o-1));
    FeatInfo{v}.CVInd(2) = str2double(stoCV(o+1:end));
    FeatInfo{v}.matname = [nam ext];

    if FeatInfo{v}.id ~= datid
        error('\n\tID Match: ID of feature extraction info does not match current data structure ID!');
    else
        fprintf('\n\tID Match: %s is ok.', nam);
    end 
end