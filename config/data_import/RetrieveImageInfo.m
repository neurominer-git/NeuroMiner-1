function [ IO, mess ] = RetrieveImageInfo(IO, datasource, mess)

if ~isfield(IO,'V') || isempty(IO.V); return; end
if ~exist('mess','var'), mess=[]; end
%Extract filenames
if iscell(IO.PP), PP = char(IO.PP); else, PP = IO.PP; end
IO.F = spm_str_manip(PP,'t'); 
switch datasource
    case {'spm','nifti'}
        if isfield(IO,'Vnan')
            IO.VV = cellcat([IO.V IO.Vnan]); 
        else
            IO.VV = cellcat(IO.V); 
        end
        L = size(IO.VV,1);
        IO.Vinfo = zeros(L,12);
        IO.Vvox = zeros(L,4);
        for i=1:size(IO.VV,1)
            IO.Vinfo(i,:) = spm_imatrix(IO.VV(i).mat);
            IO.Vvox(i,:) = diag(IO.VV(i).mat);
        end
    case 'surf'
end
