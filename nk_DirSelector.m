function D = nk_DirSelector(titlestr, D)

global SPMAVAIL

if isempty(SPMAVAIL), SPMAVAIL = logical(exist('spm_select','file')); end
if ~exist('D','var') || ~exist(D,'dir'), 
    D = {pwd}; 
elseif ischar(D)
    D = {D}; 
end
if ~exist('titlestr','var') || isempty(titlestr)
    titlestr = 'Select Directory';
end

if SPMAVAIL
    D = spm_select(1, 'dir', titlestr, D, pwd);
else
    D = uigetdir(char(D), titlestr);
end