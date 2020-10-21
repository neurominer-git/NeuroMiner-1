function PZ = nk_CheckComprFile(P)

filenotfound = fullfile(pwd,sprintf('NM_FileNotFound_%s.log',timestampstr));
filecompr = fullfile(pwd,'recompress.sh');

if iscell(P)
    PZ = cell(size(P));
    for i=1:numel(P)
        [PZ{i}] = CheckComprFile(P{i},filenotfound, filecompr);
    end
elseif size(P,1)>1
    PZ = []; MZ = [];
    for i=1:size(P,1)
        [PZ] = char(PZ,CheckComprFile(deblank(P(i,:)),filenotfound, filecompr));
    end
    PZ(1,:) = [];
else
    PZ = CheckComprFile(P, filenotfound, filecompr);
end

end

function PZ = CheckComprFile(P, filenotfound, filecompr)

PZ = P;
[pth,nam,ext] = fileparts(P);
comprcmd = ''; comprext = ''; targfile = filecompr;
if ~exist(P,'file')
   % if file does not exist check if a compressed version
   % exists
   PiGZ = fullfile(pth,[nam ext '.gz']);
   PiBZ = fullfile(pth,[nam ext '.bz']);
   PiZIP = fullfile(pth,[nam ext '.zip']);
   
   if exist(PiGZ,'file')
       fprintf('\nGunzip file: %s', PiGZ);
       gunzip(PiGZ)
       comprcmd = 'gzip -9'; comprext = '.gz';
       
   elseif exist(PiBZ,'file')
       % try Linux command
       fprintf('\nBunzip file: %s', PiBZ);
       system(['bzip -dv ' PiBZ]);
       comprcmd = 'bzip -9'; comprext = '.bz';
   
   elseif exist(PiZIP,'file')
       fprintf('\nUnzip file: %s', PiZIP);
       unzip(PiBZ);
       comprcmd = 'zip'; comprext = '.zip';
   end
   targfile = filenotfound;
elseif strcmp(ext,'.nii.gz') || strcmp(ext,'.img.gz')
    gunzip(P);
    comprcmd = 'gzip -9'; comprext = '.gz';
        
elseif strcmp(ext,'.nii.bz') || strcmp(ext,'.img.bz')
    % Works only for linux
    system(['bzip -d ' P]);
    comprcmd = 'bzip -9'; comprext = '.bz';
    
elseif strcmp(ext,'.nii.zip') || strcmp(ext,'.img.zip')
    unzip(P);
    comprcmd = 'zip'; comprext = '.zip';
end

if ~strcmp(comprcmd,'')
    ext = regexprep(ext,comprext,'');
    PZ = fullfile(pth,[nam ext]);    
    PZ = regexprep('/',PZ,'//');
    cmd = sprintf('\n%s %s', comprcmd, PZ); write2file(targfile,cmd);
end

end

function write2file(fil, txt)

fid = fopen(fil,'a');
fprintf(fid,txt);
fclose(fid);

end