function sVimg = nk_SmoothImage(Vdat, Vm, indvol, SmoothKern)

% if exist('badcoords','var') && ~isempty(badcoords), 
%     indvol = indvol(~badcoords); 
% end

Vol             = Vm;
Vol.dim         = [Vm.dim(1:3)];
Vol.dt          = [64,0];
Vimg            = zeros(Vm.dim(1:3));
sVimg           = Vimg;
Vimg(indvol)    = Vdat;
sVimg           = spm_smooth(Vimg,sVimg,SmoothKern);
sVimg           = reshape(sVimg,prod(size(sVimg)),1);
end