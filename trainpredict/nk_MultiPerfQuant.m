function Perf = nk_MultiPerfQuant(expected, predicted, modus)

nsubj = numel(expected);
ngroups = unique(expected(~isnan(expected)));
inanp = isnan(predicted);
inane = isnan(expected);
nG = numel(ngroups);

switch modus
    case 0
        errs = predicted ~= expected;
        Perf= ( 1 - sum(errs)/ nsubj) * 100;
    case 1
        Perfi = zeros(nG,1); 
        for i = 1:nG
            expi=-1*ones(size(expected)); predi=-1*ones(size(expected));
            indp = predicted == ngroups(i); inde = expected == ngroups(i);
            expi(inde)=1; predi(indp)=1;
            expi(inane)=NaN; predi(inanp)=NaN;
            Perfi(i) = BAC(expi,predi);
        end
        Perf = mean(Perfi);
end

end