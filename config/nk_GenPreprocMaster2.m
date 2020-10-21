function [featmat, emptfl, nvar] = nk_GenPreprocMaster2(datid, cv, mastermat, tdir, basename, savefl, varind, varstr, concatfl, UIopt)
% 
% This function assembles a bunch of paths to preprocessed data files (each file
% contains the data of a CV2 data partition) into a single master file.
% It allows to combine the paths to files of different variates to be used
% by the subsequent classifier / predictor training routines (see e.g.
% nk_Grid)
%
% This function belongs to the NeuroMiner software suite
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Nikolaos Koutsouleris, 05/2015

global SAV 
featmat = [];
saveflx = false; if exist('savefl','var') && ~isempty(savefl), if savefl, saveflx = true; end; end
if ~exist('tdir','var') || isempty(tdir) || ~exist(tdir,'dir'), tdir = pwd; end
emptfl = false;
if ~exist('concatfl','var') || isempty(concatfl), concatfl = false; end
if (~exist('mastermat','var') || ~exist(mastermat,'file')), 
    masterflx = false;
else
    masterflx = true;
end

if ~exist('UIopt','var') || isempty(UIopt),
    UIopt.title           = 'Generate PreprocMaster using ...';
    UIopt.filtermaster    = 'PreprocMaster.*ID.*\.mat';
    UIopt.promptmaster    = 'PreprocMaster';
    UIopt.filtermanual    = 'PreprocData.*ID.*\.mat';
    UIopt.promptmanual    = 'PreprocData';
end
fl = 1;
if ~masterflx && ~saveflx
    fl = nk_input(UIopt.title,0,'mq', ...
                'Master file|Manual file(s) selection',[1,2],1);
    switch fl
        case 0
            return
        case 1, 
            mastermat = nk_FileSelector(1, 'matrix', ['Select ' UIopt.promptmaster], UIopt.filtermaster,[],tdir);
            if isempty(mastermat) || ~exist(deblank(mastermat),'file'), return; end
    end
elseif ~masterflx 
    fl = 2;
end

switch fl
    case 1
        master=load(mastermat);
        if master.datid ~= datid
             error(['\nID Match: ID of ' UIopt.promptmaster ' file does not match ' ... 
                    ' current data structure ID!']);
        end
        featmat = master.featmat;
        nvar = master.nvar;
        
    case 2
        
        [ix, jx] = size(cv.TrainInd);
        
        % Ask for more than one variates
        nvar = 1; if exist('varind','var'), nvar = length(varind); end
        nfiles = zeros(nvar,1);
        FeatInfo = cell(nvar,1); featmat = cell(nvar,1);
        promptstr = ['Select ' UIopt.promptmanual ];
        PP = cellstr(nk_FileSelector(Inf, 'matrix', promptstr, UIopt.filtermanual,[],tdir));
        if strcmp(PP,''), emptfl = true; return; end

        for i=1:nvar % Loop through the variates
            
            if ~exist('varstr','var') || isempty(varstr)
                varstri = sprintf('var%g',varind(i));
            else
                varstri = sprintf('%s',varstr);
            end
            
            indi = find(~cellfun(@isempty,regexp(PP,varstri)));
            if isempty(indi)
                 indi = find(~cellfun(@isempty,regexp(PP,'Meta')));
            end
            nfiles(i) = numel(indi);
            FeatInfo{i} = cell(nfiles(i),1); sortarr = [];
            featmat{i} = PP(indi);
            
            for v=1:length(featmat{i})

                [pth,nam,ext] = fileparts(deblank(featmat{i}{v}));
                
                % Read ID string
                FeatInfo{i}{v}.id = nk_ExtractID(nam, ext);
                
                % Check if ID string had correct format
                if isempty(FeatInfo{i}{v}.id)
                    fprintf('\n\tNo ID found in filename %s', nam);
                    fprintf('\n\tOpening mat file to find ID.');
                    load(deblank(featmat{i}{v}));
                    try 
                        FeatInfo{i}{v}.id = id;
                    catch
                        error('\nNo ID found. Abort!')
                    end
                end
                
                    % Check if ID matches current ID of data structure
                if ~strcmp(FeatInfo{i}{v}.id,datid)
                    error('\nID Mismatch: ID of feature extraction info (%s) does not match current NM workspace ID (%s)!',FeatInfo{i}{v}.id,datid);
                else
                    fprintf('\n\tID Match: %s is ok.', nam);
                end 
                
                % Extract ID number
                FeatInfo{i}{v}.id = char(FeatInfo{i}{v}.id);
                FeatInfo{i}{v}.id = FeatInfo{i}{v}.id(3:end-1);
                FeatInfo{i}{v}.matname = fullfile(pth,[nam ext]);
                
                % Load CV2 partition info
                load(FeatInfo{i}{v}.matname,'outfold','outperm')
                sortarr = [sortarr; outperm outfold];
            end
            
            % Sort files first according to their dimensionality, their
            % outer and their inner CV index
            [S, I] = sortrows(sortarr);
            FeatInfo{i} = FeatInfo{i}(I);          
            featmat{i} = cell(ix, jx);
            for v=1:size(FeatInfo{i},1)
                featmat{i}(S(v,1), S(v,2)) = {FeatInfo{i}{v}.matname};
            end
            
        end
        
        % Concatenate file paths if more than one variate is available
        if nvar > 1 && concatfl
            % Concatenate paths
            fprintf('\nConcatenating paths of modality data.')
            F = cell(ix, jx);
            for j= 1:ix
                for k=1:jx
                    for l=1:nvar
                         F{j,k} = char(F{j,k},featmat{l}{j,k});
                    end
                    F{j,k}(1,:)=[];
                end
            end 
            featmat = []; featmat{1} = F;
        end
        
        % Check the saving options
        nFsel = ~cellfun(@isempty,featmat{1}); 
        if (~exist('savefl','var') || isempty(savefl)) && sum(nFsel(:))>1
            savefl = nk_input(['Save as ' UIopt.promptmaster ' file?' ],0,'yes|no',[1,0],1);
        end
        
        % if saving is activated, go on to save the master path file
        if savefl
            if ~exist('tdir','var') || ~exist(tdir,'dir')
                tdir = nk_DirSelector('Select a target directory for the master file(s)');
            end
            if ~exist('basename','var')
               basename = [nk_input('Specify a basename for the master files',0,'s', ...
                   SAV.matname) '_']; 
            end
            savpth = fullfile(tdir,[basename UIopt.promptmaster '_ID' datid ext]); 
            save(savpth,'featmat','datid', 'nvar');
        end     

end