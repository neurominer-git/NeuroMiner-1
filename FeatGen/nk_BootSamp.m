function BootSamp = nk_BootSamp(labels,nboot,perc)

m = length(labels);
minl=min(labels);
maxl=max(labels);
BootSamp = zeros(m,nboot);

for i=1:nboot
    for j=minl:maxl
        ind = find(labels==j);
        lind = length(ind);
        pl = round(lind*perc/100);
        pind = randsample(ind,pl,false);
        
        %BootSamp(:,i) = x(floor(size(x,1)*rand(size(x,1),1))+1,:); % bootstrap observations
        BootSamp(ind,i) = randsample(pind,lind,true);
    end
    %uq(i) = length(unique(BootSamp(:,i)))*100/m;
end

%fprintf('\nUniques: %g',mean(uq));

return