function sV = nk_VisSmooth(V, fwhm, varind, brainmask, badcoords, dimvecx)

sV = zeros(size(V));

for i=1:numel(varind)
    
    if numel(varind)>1
        fprintf('\n\nSmoothing weight vector data for Variate #%g',i)
    end
    % Get pointer2variate vector
    if exist('dimsizes','var') 
        varsuff = sprintf('_var%g',varind(i));
        brainmaski = brainmask{i};
        if ~isempty(badcoords), 
            if iscell(badcoords),
                badcoordsi = badcoords{i}; 
            else
                badcoordsi = badcoords;
            end
        else
            badcoordsi=[]; 
        end    
    end
    dimvec = dimvecx(i)+1:dimvecx(i+1);
    
    % Reshape vector into 3D image
    [Vi, indvol, VOX] = nk_WriteVol(V(dimvec),[],2,brainmaski,badcoordsi);
    
    % Smooth 3D image
    Q = zeros(size(Vi));
    fwhm=fwhm./VOX;
    spm_smooth(Vi,Q,fwhm);
    vQ = Q(indvol);
    sV(dimvec) = vQ;
    
end