function DeleteTempVol(V, flag)
if exist(V.fname,'file') && flag
    [p,n] = fileparts(V.fname);
     delete(V.fname); 
     hdrfile = fullfile(p,[n '.hdr']);
     if exist(hdrfile,'file'), delete(hdrfile); end
end