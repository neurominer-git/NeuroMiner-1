function [featmat, emptfl, skiparr, missarr] = nk_GenCVdataMaster2(datid, cv, mastermat, tdir, basename, savefl, varind, varstr, concatfl, UIopt)

global SAV
featmat = []; skiparr=[]; missarr=[];

saveflx = false; if exist('savefl','var') && ~isempty(savefl), if savefl, saveflx = true; end; end
if ~exist('tdir','var') || isempty(tdir) || ~exist(tdir,'dir'), tdir = pwd; end
emptfl = false;
if (~exist('mastermat','var') || ~exist(mastermat,'file')), 
    masterflx = false;
else
    masterflx = true;
end
if ~exist('concatfl','var') || isempty(concatfl), concatfl = false; end

if ~exist('UIopt','var') || isempty(UIopt),
    UIopt.title           = 'Generate CVdatamat master file using ...';
    UIopt.filtermaster    = 'CVdataMaster.*ID.*\.mat';
    UIopt.promptmaster    = 'CVdataMaster';
    UIopt.filtermanual    = 'CVdatamat.*ID.*\.mat';
    UIopt.promptmanual    = 'CVdatamat';
end

if ~masterflx && ~saveflx
    fl = nk_input(UIopt.title,0,'mq', ...
                'Master file|Manual file(s) selection',[1,2],1);
    switch fl 
        case 0
            return
        case 1
            mastermat = nk_FileSelector(1, 'matrix', ['Select ' UIopt.promptmaster], UIopt.filtermaster,[],tdir);
            if isempty(mastermat) || ~exist(deblank(mastermat),'file'), return; end
    end
elseif ~masterflx
    fl = 2;
else
    fl = 1;
end

switch fl
    
    case 1
 
        master=load(mastermat);
        if master.datid ~= datid
             error(['\nID Match: ID of ' UIopt.promptmaster ' file does not match ' ... 
                    ' current data structure ID!']);
        end
        featmat = master.featmat;

    case 2
        
        [ix, jx] = size(cv.TrainInd);
        % Ask for more than one variates
        nvar = 1; if exist('varind','var'), nvar = numel(varind); end
        nfiles = zeros(nvar,1);
        FeatInfo = cell(nvar,1); featmat = cell(nvar,1);
        promptstr = ['Select ' UIopt.promptmanual ];
        PP = cellstr(nk_FileSelector(Inf, 'matrix', promptstr, UIopt.filtermanual,[],tdir));
        if strcmp(PP{1},''), emptfl = true; return; end        
        
        for i=1:nvar % Loop through the variates
            
            if ~exist('varstr','var') || isempty(varstr)
                varstri = sprintf('var%g',varind(i));
            else
                varstri = varstr;
            end
            
            indi = find(~cellfun(@isempty,regexp(PP,varstri)));
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
                    load(deblank(featmat{i}{v}),'id');
                    try 
                        FeatInfo{i}{v}.id = id;
                    catch
                        [~, filnamex] = fileparts( deblank(featmat{i}{v}) );
                        II = strfind(filnamex,'_ID');
                        if ~isempty(II)
                            FeatInfo{i}{v}.id = filnamex(II+3:numel(filnamex));
                        else
                            error(sprintf('\nNo ID found. Abort!'))
                        end
                    end
                end
                
                    % Check if ID matches current ID of data structure
                if ~strcmp(FeatInfo{i}{v}.id,datid)
                    warning('\nID Mismatch: ID of feature extraction info (%s) does not match current NM workspace ID (%s)!',FeatInfo{i}{v}.id,datid);
                else
                    fprintf('\n\tID Match: %s is ok.', nam);
                end 
                
                % Extract ID number
                FeatInfo{i}{v}.id = char(FeatInfo{i}{v}.id);
                FeatInfo{i}{v}.id = FeatInfo{i}{v}.id(3:end-1);
                FeatInfo{i}{v}.matname = fullfile(pth,[nam ext]);
                
                 % Load CV2 partition info
                try
                    % Load CV2 partition info
                    load(FeatInfo{i}{v}.matname,'ofold','operm')
                    sortarr = [sortarr; ofold operm];
                catch
                    warning(['Cannot open ' nam '. Skip']);
                    skiparr = char(skiparr, deblank(featmat{i}{v}));
                end
            
            end
            if ~isempty(skiparr)
                fprintf('\n\nCould not open the following files:')
                for u=1:size(skiparr,1)
                    fprintf('\n%s',skiparr(u,:))
                end
                abortfl = nk_input('What to do?',0,'m','Abort|Erase damaged files and abort|Erase damaged files',[1,2,3]);
                switch abortfl
                   case 1
                       return
                   case {2,3}
                        for u=1:size(skiparr,1)
                            delete(deblank(skiparr(u,:)))
                        end 
                        if abortfl == 2, return; end
                end
            end
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
