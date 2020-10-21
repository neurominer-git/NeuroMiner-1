function [V, flag] = WriteTempVol(V,Vvol, flag)

if ~exist(V.fname,'file')
    [~,n,e] = fileparts(V.fname);
    tmpfile = fullfile(pwd,[n e]);
    
    V.fname = tmpfile;
    spm_write_vol(V, Vvol)
    flag = true;
end