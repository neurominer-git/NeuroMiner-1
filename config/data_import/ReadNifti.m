function IO = ReadNifti(IO)

% Scale images
if isfield(IO,'globscale') && IO.globscale>1
    g = IO.g;
    GM = IO.globnorm;
else
    g = ones(size(IO.PP,1),1); GM = 1;
    cprintf('blue','\nNo global normalization')
end

gSF = GM./g; 
V = IO.V;
cnt = 1; nV = numel(V); 

for j = 1 : nV
    for i = 1:numel(V{j})
        V{j}(i).pinfo(1:2,:) =  V{j}(i).pinfo(1:2,:) * gSF(cnt);
        cnt = cnt+1;
    end
end

if isfield(IO,'nangroup') 
    if IO.nangroup && IO.nan_subjects > 0
        V{nV+1} = IO.Vnan;
        for i=1:IO.nan_subjects
            V{nV+1}(i).pinfo(1:2,:)  =  V{nV+1}(i).pinfo(1:2,:) * gSF(cnt);
            cnt = cnt + 1;
        end
    end
end
%% Print image import info
if iscell(IO.PP), IO.PP = char(IO.PP); end
N = size(IO.PP,1);
F = cellstr([ spm_str_manip(IO.PP,'c') ...
                        repmat(' | Global: ',N,1) num2str(g,'%4.1f') ...
                        repmat(' | Scaling: ',N,1) num2str(gSF,'%g') ...
                        repmat(' | Vox dim: ',N,1) num2str(IO.Vvox,'%1.2f ') ...
                        repmat(' | Origin: ',N,1) num2str(IO.Vinfo(:,1:3), '%1.1f' )]);
filereport = fullfile(pwd,'ReadNifti_Report.txt');
fid=fopen(filereport,'a');
fprintf(fid,'NeuroMiner Report: Import NiFTI/Analyze files ');

if IO.oocvflag,
    fprintf(fid,'\nCaller: Independent test data import');
else
    fprintf(fid,'\nCaller: Discovery data import');
end

fprintf(fid,'\nDate: %s', datestr(now));
fprintf(fid,'\nOS: %s',computer);
fprintf(fid,'\n%s', repmat('=',1,size(char(F),2)));
fprintf(fid,'\n%s', F{:});
fprintf(fid,'\n\n');
fclose(fid);
clc
cprintf('black*','See %s for import information \n', filereport);

if IO.nangroup && IO.nan_subjects > 0, 
    n_subjects = [IO.n_subjects IO.nan_subjects];
else
    n_subjects = IO.n_subjects;
end

tmpflg = false;
if ~isfield(IO,'Vmvol'), IO.Vmvol = spm_read_vols(IO.Vm); end

% Read images for weighting / subspace extraction
if isfield(IO,'Vw') && ~isempty(IO.Vw)
    
    [IO.Vm, tmpflg] = WriteTempVol(IO.Vm, IO.Vmvol, tmpflg);
    IO.Yw = nk_ReturnSubSpaces(IO.Vw, IO.Vm, 1, numel(IO.Vw) , IO.Thresh, filereport);
end

% Read images
[IO.Vm, tmpflg] = WriteTempVol(IO.Vm, IO.Vmvol, tmpflg);
IO.Y = nk_ReturnSubSpaces(V, IO.Vm, numel(V), n_subjects, IO.Thresh, filereport);

DeleteTempVol(IO.Vm, tmpflg) 




