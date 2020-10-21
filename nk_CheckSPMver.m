function spmver = nk_CheckSPMver

spmver=0;
if ~isempty(which('spm')) && (strcmp(spm('ver'),'SPM5') || strcmp(spm('ver'),'SPM8')|| strcmp(spm('ver'),'SPM12')),  
    spmver = 1; 
end
