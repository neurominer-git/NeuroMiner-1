function [featmat,emptfl] = nk_GenCVresultsMaster(datid, mastermat, tdir, basename, savefl)

global SAV

if exist('savefl','var') && savefl, saveflx = true; else saveflx = false; end
if ~exist('tdir','var') || isempty(tdir) || ~exist(tdir,'dir'), tdir = pwd; end
emptfl = false;
if (~exist('mastermat','var') || ~exist(mastermat,'file')), masterflx = false; else masterflx = true; end
fl = 1;
if  ~masterflx && ~saveflx
    fl = nk_input('Generate CVresultsMaster using ...',0,'mq', ...
    'CVresultsMaster file|Manual CVresults file(s) selection',[1,2],1);
    switch fl
        case 0
            return
        case 1
            mastermat = nk_FileSelector(1, 'matrix', 'Select CVresults Master', 'CVresultsMaster.*ID.*\.mat',[],tdir);
            if isempty(mastermat) || ~exist(deblank(mastermat),'file'), return; end
    end
elseif ~masterflx
    fl = 2;
end

switch fl 

    case 1
        master=load(mastermat);
        if master.datid ~= datid
             error(['\nID Match: ID of CVresults master file does not match ' ... 
                    ' current data structure ID!']);
        end
        featmat = master.featmat;
    case 2
        gdmat = nk_FileSelector(Inf, 'matrix', 'Select CVresults', '.*CVresults.*_ID.*\.mat',[],tdir);
        if strcmp(gdmat,''), emptfl = true; return; end 
        inds = zeros(size(gdmat,1),1);
    
        for i=1:size(gdmat,1)
            [~,nam,ext] = fileparts(deblank(gdmat(i,:)));
            inds(i) = str2double(char(regexp(nam, '_var(\d+)_$*', 'tokens', 'once')));
            id = nk_ExtractID(nam, ext);
            if ~strcmp(id,datid)
                warning('ID Match: ID of feature extraction info does not match current data structure ID! Skipping!');
                continue;
            else
                fprintf('\n\tID Match: %s is ok.', nam);
            end
        end
        [~,inds] = sort(inds);
        featmat = cellstr(gdmat(inds,:));
        if numel(featmat)>1
            nFsel = ~cellfun(@isempty,featmat{1}); 
        end
        if ~exist('savefl','var') || isempty(savefl) && sum(nFsel(:))>1
            savefl = nk_input('Save as CVresultsMaster file?',0,'yes|no',[1,0],1);
        end
       
        if savefl
            
            if ~exist('tdir','var') || ~exist(tdir,'dir')
               tdir = nk_DirSelector('Select a target directory for the master file(s)');
            end

            if ~exist('basename','var')
               basename = [nk_input('Specify a basename for the master files',0,'s', ...
                   SAV.matname) '_']; 
            end
            
            savpth = fullfile(tdir,[basename 'CVresultsMaster_ID' datid ext]); 
            save(savpth,'featmat','datid');
        end   
end

return
