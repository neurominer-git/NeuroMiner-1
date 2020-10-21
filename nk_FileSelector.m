function [imgfile, imgvol, mess] = nk_FileSelector(nimg, imgtype, titlestr, fltstr, imgfile, startdir, mess)
% =========================================================================================
% [imgfile, imgvol, mess] = nk_FileSelector(nimg, imgtype, titlestr, fltstr, imgfile, mess)
% =========================================================================================
% NeuroMiner file selection function using either the SPM file selection
% tool or the MATLAB builtin uigetfile, depending on the availability of
% SPM in NM. 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Nikolaos Koutsouleris, 03/2017

global SPMAVAIL

if ~exist('mess','var'), mess=[]; end
if ~exist('startdir','var') || ~exist(startdir,'dir'), startdir = pwd; end
imgvol = [];
if isempty(SPMAVAIL), SPMAVAIL = logical(exist('spm_select','file')); end
if exist('imgfile','var') && ~isempty(imgfile) && ischar(imgfile)
    if size(imgfile,1)>1
        imgfile = char(regexprep(cellstr(imgfile),',1',''));
    else
        imgfile = regexprep(imgfile,',1','');
    end
else
    imgfile = [];
end

if SPMAVAIL
    if ~isempty(imgfile) && ~iscell(imgfile),
        imgfile = cellstr(imgfile);
    end
    imgfile = spm_select(nimg, fltstr, titlestr, imgfile, startdir, fltstr);
else
    fltstr = regexprep(fltstr,'\\.',''); fltstr = regexprep(fltstr,'\.',''); 
    if ~strcmp(fltstr(1),'*'), fltstr = ['*' fltstr]; end
    [imgfile, imgpath] = uigetfile(fltstr,titlestr, startdir, 'MultiSelect', 'on');
    if iscell(imgfile)
        imgfile = cellfun(@fullfile, repmat({imgpath},1,numel(imgfile)), imgfile,'UniformOutput',false);
        imgfile = char(imgfile');
    else
        imgfile = fullfile(imgpath,imgfile);
    end
end

if isempty(imgfile) && (~ischar(imgfile) || iscell(imgfile)), return; end

switch imgtype
    case {'nifti','spm','surf'}
        switch imgtype
            case {'nifti','spm'}
                imgfile = nk_CheckComprFile(imgfile);
        end
        [imgvol, imgfile, mess] = ReadVol(imgfile, imgtype, mess);
end  

function [imgvol, imgfile, mess] = ReadVol(imgfile, imgtype, mess)

if isempty(imgtype), imgvol=[]; return; end
ind = true(size(imgfile,1),1);
imgvol = cell(size(imgfile,1),1);
fprintf('\n');
for i=1:size(imgfile,1)
    fprintf('Reading %s\r',imgfile(i,:));
    try
        iimgfile = deblank(imgfile(i,:));
        switch imgtype
            case {'nifti','spm'}
                timgvol = spm_vol(iimgfile); 
            case 'surf'     
                timgvol = SurfaceReader(iimgfile);
        end
        clear timgvol
    catch
        ind(i) = false;
        mess = GenerateMessageEntry(mess, sprintf('ERROR: File %s could not be opened and was removed from the file list! ',imgfile(i,:)));
        if ~exist('fid','var')
            fid = fopen(fullfile(pwd,sprintf('NM_ReadError_%s.log',timestampstr)),'a');
            fprintf(fid,'The following image files could not be opened:');  
        else
            fprintf(fid,'\n%s',imgfile(i,:));
        end
    end
end
if exist('fid','var'), fclose(fid); end
imgfile = imgfile(ind,:);
if ~isempty(imgfile)
    switch imgtype
        case {'nifti','spm'}
            imgvol = spm_vol(imgfile); 
        case 'surf'
            if size(imgfile,1)==1,
                iimgfile = deblank(imgfile(1,:));
                imgvol = SurfaceReader(iimgfile);
            else
                imgvol = cell(size(imgfile,1),1);
                for i=1:size(imgfile,1)
                    iimgfile = deblank(imgfile(i,:));
                    [~,~,e] = fileparts(iimgfile);
                    imgvol{i} = SurfaceReader(iimgfile);
                end
            end
    end
end