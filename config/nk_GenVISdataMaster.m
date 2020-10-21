function [featmat, emptfl, skiparr, missarr] = nk_GenVISdataMaster(datid, cv, mastermat, tdir, basename, savefl)

global SAV
featmat = []; skiparr=[]; missarr=[];
saveflx = false; if exist('savefl','var') && savefl, saveflx = true; end
if (~exist('mastermat','var') || ~exist(mastermat,'file')), 
    masterflx = false;
else
    masterflx = true;
end
emptfl = false;
if ~masterflx && ~saveflx
    fl = nk_input('Generate VISdataMaster using ...',0,'mq', ...
    'VISdataMaster file|Manual VISdatamat file(s) selection',[1,2],1);
    if ~fl
        return
    elseif fl == 1
        mastermat = nk_FileSelector(1, 'matrix', 'Select VISdatamat Master', 'VISdataMaster.*ID.*\.mat');
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
             error(['\nID Match: ID of VISdatamat master file does not match ' ... 
                    ' current data structure ID!']);
        end
        featmat = master.featmat;

    case 2
        
        [ix, jx] = size(cv.TrainInd);
        gdmat = nk_FileSelector(Inf, 'matrix', 'Select VISdatamats', '*VISdatamat.*ID.*\.mat');
        if strcmp(gdmat,''), emptfl = true;return; end
        sortarr=[]; 
        
        for i=1:size(gdmat,1)
            
            [~,nam,ext] = fileparts(deblank(gdmat(i,:)));
            
            id = char(regexp([nam ext],'ID.*\.','match'));
            id = id(3:end-1);
            
            if id ~= datid
                warning('ID Match: ID of feature extraction info does not match current data structure ID! Skipping!');
                continue;
            else
                fprintf('\n\tID Match: %s is ok.', nam);
            end
            
            % Load CV2 partition info
            try
                load(deblank(gdmat(i,:)),'ofold','operm')
                sortarr = [sortarr; ofold operm];
            catch
                warning(['Cannot open ' nam '. Skip']);
                skiparr = char(skiparr, deblank(gdmat(i,:)));
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
        gdmat = gdmat(I,:);
        featmat = cell(ix, jx);
        for i=1:size(gdmat,1)
            featmat(S(i,2), S(i,1)) = {deblank(gdmat(i,:))};
        end
        missarr = cellfun('isempty',featmat);
       
        if ~exist('savefl','var') || isempty(savefl)
            savefl = nk_input('Save as VISdata-Master file?',0,'yes|no',[1,0],1);
        end
        
        if savefl
            
            if ~exist('tdir','var') || ~exist(tdir,'dir')
                 tdir = nk_DirSelector('Select a target directory for the master file(s)');
            end
            cd(tdir);

            if ~exist('basename','var')
               basename = [nk_input('Specify a basename for the master files',0,'s', ...
                   SAV.matname) '_']; 
            end
            savpth = fullfile(tdir,[basename 'VISdataMaster_ID' datid ext]); 
            save(savpth,'featmat','datid');
        end
end


return
