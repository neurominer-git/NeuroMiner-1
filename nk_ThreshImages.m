function nk_ThreshImages

P = spm_select(Inf,'image','Select images');
V = spm_vol(P);

threshval = spm_input('Specify percentile threshold for images ',0,'e');

for i=1:numel(V)
    [p,n,e] = fileparts(V(i).fname);
    Vx = spm_read_vols(V(i));
    ind = Vx(:) > 0; 
    Thresh = prctile(Vx(ind), threshval);
    VxP = zeros(size(Vx));
    indT = Vx >= Thresh;
    VxP(indT) = Vx(indT);
    Vw = V(i);
    Vw.fname = fullfile(p,[n '_prctile-gr' num2str(threshval) e]);
    spm_write_vol(Vw,VxP);
end
