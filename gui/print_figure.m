function print_figure(handles, objnames, savename)

tmpfig = figure;
%tmpfig = figure('visible','off');
n_obj = numel(objnames);
h = [];
for i= 1:n_obj
    hi = handles.(objnames{i});
    h =[hi h];  
end
newax = copyobj(h,tmpfig); 
set(newax, 'units', 'normalized', 'position', [0.13 0.11 0.775 0.815]);
%I = getframe(tmpfig);
%imwrite(I.cdata, savename);
hgsave(tmpfig,savename)
close(tmpfig)