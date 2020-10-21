function FEAT = nk_PerfCommonFeatExtraction(FEAT, param)
global PREPROC

[iy, jy] = size(FEAT);
n = iy*jy;
[d, nclass] = size(FEAT{1,1}.dat{1}.vol);
switch param.mode
    case {1,2,3}
        Dx = false(d, nclass);
        param.mode = 1;
        param.thresh = 75;
    case 4
        Dx = zeros(d, nclass);
    case 5
        try 
            percvec = PREPROC.FEATSEL.salAdap.perc_vec';
        catch
            error('Cross-CV1 threshold computation not feasible! Check percentage step vector')
        end
        Perc = zeros(nclass,1);
        nsteps = numel(percvec);
        Df = zeros(nsteps, nclass);
end
for curclass = 1:nclass

    switch param.mode
        case 1
            Df = false(d, n);
            for k = 1:n
               Df(:,k) = FEAT{k}.dat{1}.vol(:,curclass);
            end
            [P, sind] = sort(sum(Df,2) * 100 / n,'descend');
            thresh = param.thresh;
        case 2
            Df = zeros(d, n);
            for k = 1:n
               Df(:,k) = FEAT{k}.dat{2}.vol(:,curclass);
            end
            [P, sind] = sort(mean(Df,2),'descend');
            thresh = percentile(P, param.thresh);
        case 3
            Df = zeros(d,n);
            for k = 1:n
               Df(:,k) = FEAT{k}.dat{1}.vol(:,curclass);
            end
            thresh = mean(sum(Df));
        case 4
            fprintf('\nComputing feature selection probability vector #%g.',curclass)
            Df = zeros(d,1);
             for k = 1:n
               Df = Df + FEAT{k}.dat{1}.vol(:,curclass);
             end
            Dx(:,curclass) = Df / n;
        case 5

            for k = 1:n
                if size(FEAT{k}.YPerf{curclass},2) > 1,
                    Yperf = FEAT{k}.YPerf{curclass}';
                else
                    Yperf = FEAT{k}.YPerf{curclass};
                end
                Df(:,curclass) = Df(:,curclass) + Yperf;
            end
            Df(:,curclass) = Df(:,curclass)/n;    
            [opt, indopt] = max(Df(:,curclass));
            Perc(curclass) = percvec(indopt);
            fprintf('\nFound optimum = %g for predictor %g at %g%% percentile', opt, curclass, Perc(curclass))
            
    end
    switch param.mode
        case {1,2,3}
            ind = P > thresh;
            fprintf('\nSelected %g features > consistency threshold %1.2f%% for dichotomization #%g', sum(ind), thresh, curclass)
            Dx(sind(ind),curclass) = true;
    end
end

for curclass = 1 : nclass
    for k = 1:n
        switch param.mode
            case {1,2,3}
                FEAT{k}.dat{1}.vol(:,curclass) = FEAT{k}.dat{1}.vol(:,curclass) & Dx(:,curclass);
            case 4
                FEAT{k}.dat{2}.vol(:,curclass) = Dx(:,curclass);
            case 5
                thresh = percentile(FEAT{k}.dat{2}.vol(:,curclass),Perc(curclass));
                FEAT{k}.dat{1}.vol(:,curclass) = FEAT{k}.dat{2}.vol(:,curclass) > thresh;
                
        end
    end
end

end