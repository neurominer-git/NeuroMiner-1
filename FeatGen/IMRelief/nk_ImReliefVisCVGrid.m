function Res = nk_ImReliefVisCVGrid(TempDir, CV, FEATSEL)

if ~exist('TempDir','var') || isempty(TempDir)
    RootDir = spm_select(1,'dir','Select IMRelief root directory');
    if isempty(RootDir), reutrn; end
end

[niy, njy] = size(CV.TrainInd);
[nix, njx] = size(CV.cvin{1,1}.TrainInd);
switch FEATSEL.BINMOD
    case 1
        nclass = numel(CV.class{1,1});
    case 0
        nclass = 1;
end
Res.mMxCl = cell(nclass,1);
Res.sMxCl = cell(nclass,1);
Res.bestlambda = zeros(1,nclass);
Res.bestsigma  = zeros(1,nclass);
Res.mThreshCl = cell(nclass,1);
Res.bestthresh = zeros(1,nclass);
Res.mCritCl = cell(nclass,1);
Res.bestcrit = zeros(1,nclass);

ho = findobj('Name','IMRelief grid parameter evaluation');
if isempty(ho), 
    h = figure('Name','IMRelief grid parameter evaluation','NumberTitle','off'); 
else
    h = ho;
end

algopref = GetSalModeStr(FEATSEL);

for curclass = 1:nclass
    
    if iscell(FEATSEL.imrelief)
        if numel(FEATSEL.imrelief)<nclass
            imrelief = FEATSEL.imrelief{end};
        else
            imrelief = FEATSEL.imrelief{curclass};
        end
    else
        imrelief = FEATSEL.imrelief;
    end
    
    nsigma = numel(imrelief.sigma);
    nlambda = numel(imrelief.lambda);
    MxCl = zeros(nsigma, nlambda, 1);
    ThreshCl = zeros(nsigma, nlambda,1);
    CritCl = zeros(nsigma, nlambda,1);
    Res.mMxCl{curclass} = zeros(nsigma, nlambda);
    Res.sMxCl{curclass} = zeros(nsigma, nlambda);
    Res.mThreshCl{curclass} = zeros(nsigma, nlambda);
    Res.mCritCl{curclass} = zeros(nsigma, nlambda);
    
    ll=1;
    fprintf('\nAveraging parameter grid across CV1 partitions for classifier: %s.', ...
        CV.class{1,1}{curclass}.groupdesc)
    
    zcnt = 1;
    CosSim = zeros(nsigma, nlambda, niy*njy);
    for iy = 1:niy

        for jy = 1:njy
            
            TempDir = fullfile(RootDir, ['IMRelief_oCV' num2str(iy) '.' num2str(jy)]);
            if ~exist(TempDir,'dir'), continue; end
            
            iP = allcomb(1:nix,1:njx);
            niP = length(iP); 
            
            C=zeros(nsigma,nlambda);
            cnt=0;
            fprintf('\nComputing pairwise similarity matrix of weight vectors across CV1:')
            for ix = 1:niP-1
                TempFile = fullfile(TempDir, ...
                                ['IMRelief_Eval-' algopref '_CVData_iCV' num2str(iP(ix,1)) '.' num2str(iP(ix,2)) ...
                                '_cl' num2str(curclass) '.mat']);
                if exist(TempFile,'file')
                    D1 = load(TempFile,'DScore'); 
                else
                    break
                end
                fprintf('\n')
                for jx = ix:niP
                    TempFile = fullfile(TempDir, ...
                                ['IMRelief_Eval-' algopref '_CVData_iCV' num2str(iP(jx,1)) '.' num2str(iP(jx,2)) ...
                                '_cl' num2str(curclass) '.mat']);
                    if exist(TempFile,'file')
                        D2 = load(TempFile,'DScore');
                    else
                        break
                    end
                    tSim = compute_weightvecsim(D1.DScore,D2.DScore);
                    C = C + tSim;
                    fprintf('.')
                    cnt=cnt+1;
                end
            end
            CosSim(:,:,zcnt) = C ./ cnt;
            for ix = 1:nix

                for jx = 1:njx
                    
                    flg = false;
                    
                    TempFile = fullfile(TempDir, ...
                                ['IMRelief_Eval-' algopref '_CVData_iCV' num2str(ix) '.' num2str(jx) ...
                                '_cl' num2str(curclass) '.mat']);
                    try        
                        if exist(TempFile,'file')
                            fprintf('\nLoading grid from %s', TempFile)
                            load(TempFile,'tMx', 'Thresh','Crit','dx')
                            
                            flg = true;
                            if isfield(FEATSEL,'OOCV') && ~exist('ocvfl','var')
                                load(TempFile,'Mx','Mxoocv')
                                if exist('Mx','var') && exist('Mxoocv','var')
                                    ocvfl = true;
                                    fprintf('CV-Optimization and OOCV-Optimization grid have been detected.')
                                    ocvmode = nk_input('Processing mode',0,'m', ...
                                        ['Use CV-Optimization grid|' ...
                                        'Use OOCV-Optimization grid|' ...
                                        'Use combined CV-OOCV-Optimization grid (average)|' ...
                                        'Use combined CV-OOCV-Optimization grid (product)'],1:4,1);
                                else
                                    ocvfl = false;
                                end
                            elseif ~exist('ocvfl','var')
                                ocvfl = false;
                            elseif ocvfl 
                                 load(TempFile,'Mx','Mxoocv')
                            end
                            
                        else % Try old file name
                            TempFile = fullfile(TempDir, ...
                                    ['IMRelief_CVData_iCV' num2str(ix) '.' num2str(jx) ...
                                    '_cl' num2str(curclass) '.mat']);
                            if exist(TempFile,'file')
                                fprintf('\nLoading grid from %s', TempFile)
                                load(TempFile,'tMx', 'Thresh','Crit','dx')
                                flg = true;
                            end
                        end
                        
                        if flg, 
                            
                            if exist('dx','var')
                                ndx = numel(dx);
                                
                                if ocvfl
                                    Mx = reshape(Mx, nsigma, nlambda, ndx);
                                    Mxoocv = reshape(Mxoocv, nsigma, nlambda, ndx);
                                end
                                
                                if ll == 1;
                                    MxCl = zeros(nsigma, nlambda, 1, ndx);
                                    ThreshCl = zeros(nsigma, nlambda, 1, ndx);
                                    CritCl = zeros(nsigma, nlambda, 1, ndx);
                                    Res.mMxCl{curclass} = zeros(nsigma, nlambda, ndx);
                                    Res.sMxCl{curclass} = zeros(nsigma, nlambda, ndx);
                                    Res.mThreshCl{curclass} = zeros(nsigma, nlambda, ndx);
                                    Res.mCritCl{curclass} = zeros(nsigma, nlambda, ndx);
                                end
                                for d = 1: ndx
                                    if ocvfl 
                                        switch ocvmode
                                            case 1
                                                MxCl(:, :, ll, d) = Mx(:,:,d);
                                            case 2
                                                MxCl(:, :, ll, d) = Mxoocv(:,:,d);
                                            case 3
                                                MxCl(:, :, ll, d) = (Mx(:,:,d) + Mxoocv(:,:,d)) / 2;
                                            case 4
                                                MxCl(:, :, ll, d) = (Mx(:,:,d) .* Mxoocv(:,:,d)) / 100;
                                        end
                                    else
                                        MxCl(:, :, ll, d) = tMx(:,:,d);
                                    end
                                    if exist('Thresh','var'), ThreshCl(:,:,ll, d) = Thresh; clear Thresh; end
                                    if exist('Crit','var'), CritCl(:, :, ll, d) = Crit; clear Crit; end
                                end
                            else
                                ndx = 1;
                                if ocvfl 
                                    switch ocvmode
                                        case 1
                                            MxCl(:, :, ll) = Mx;
                                        case 2
                                            MxCl(:, :, ll) = Mxoocv;
                                        case 3
                                            MxCl(:, :, ll) = (Mx + Mxoocv) / 2;
                                        case 4
                                            MxCl(:, :, ll) = (Mx .* Mxoocv) / 100;
                                    end
                                else
                                    MxCl(:,:,ll) = tMx;
                                end
                                if exist('Thresh','var'), ThreshCl(:,:,ll) = Thresh; clear Thresh; end
                                if exist('Crit','var'), CritCl(:,:,ll) = Crit; clear Crit; end
                            end
                            
                            ll=ll+1;                    
                        end       
                    catch
                        continue
                    end
                end
            end
            zcnt=zcnt+1;
        end
    end
    Res.CosSim{curclass} = mean(CosSim,3);
    Res.CosSim{curclass}
    for d = 1 : ndx
        Res.mMxCl{curclass}(:,:,d) = mean(MxCl(:,:,:,d),3);
        Res.sMxCl{curclass}(:,:,d) = std(MxCl(:,:,:,d),[],3);
        Res.mThreshCl{curclass}(:,:,d) = mean(ThreshCl(:,:,:,d),3);
        Res.mCritCl{curclass}(:,:,d) = mean(CritCl(:,:,:,d),3);
        Rmd = Res.mMxCl{curclass}(:,:,d); Res.mmMxCl(d) = mean(Rmd(:));
    end
    fprintf('\nCompleted grid averaging for classifier: %s.', CV.class{1,1}{curclass}.groupdesc)
    
    fprintf('\nFound %g data partitions.',ll-1)
    
end

xpos = cell(nclass,1); ypos=cell(nclass,1);
Res.NumC1 = ll - 1;

if ndx>1
    figure;plot(dx,Res.mmMxCl,'b*-','LineWidth',2)
    title('Average grid performance across dimensions')
    xlabel('Dimensions')
    ylabel('Performance')
end
for curclass=1:nclass
    
    if iscell(FEATSEL.imrelief)
        if numel(FEATSEL.imrelief)<nclass
            imrelief = FEATSEL.imrelief{end};
        else
            imrelief = FEATSEL.imrelief{curclass};
        end
    else
        imrelief = FEATSEL.imrelief;
    end
    
    nsigma = numel(imrelief.sigma);
    nlambda = numel(imrelief.lambda);
    
    fprintf('\nBinary comparison %g',curclass)
    
    for d = 1 : ndx
        if ndx > 1, fprintf('\n\nDimensionality %g',dx(d)), end
        MxClx = Res.mMxCl{curclass}(:,:,d);
        % Filter performance grid using average function for more stable
        % parameter selection
        try
            a = fspecial('gaussian');
            MxClx = filter2(a,MxClx);
        catch
            fprintf('\nGaussian filtering not feasible')
        end
        tmxMxClx = (max(MxClx(:)) + Res.mmMxCl(d))/2;
        fprintf('\nComputed average CV1 grid across CV1 partitions:')
        try
            for z=1:nsigma,
                fprintf('\n')
                fprintf('\t%1.1f',MxClx(z,:)); 
            end
        catch
            fprintf('\nPrinting error.')
        end
        if d > 1 
            if tmxMxClx > mxMxClx, 
                mxMxClx = tmxMxClx;
                OptGrid = MxClx ;
                Opt = max(OptGrid(:));
                bestd = d;
            end
        else
            mxMxClx = tmxMxClx;
            bestd = 1;
            Opt = max(MxClx(:));
            OptGrid = MxClx;
        end
    end
    
    if ndx > 1
        fprintf('\nBest dimensionality detected at %g dimensions.', dx(bestd))
        fprintf('\nExtracting grid from this dimensionality.')
        Res.best_MxClx{curclass} = OptGrid; %Res.mMxCl{curclass}(:,:,bestd);
        Res.best_dim(curclass) = dx(bestd);
        Res.best_dim_index(curclass) = bestd;
        Res.best_param(curclass) = Opt;
    else
        Res.best_MxClx{curclass} = OptGrid;
        Res.best_param(curclass) = Opt;
    end
    
    [xpos{curclass}, ypos{curclass}] = find(Res.best_MxClx{curclass} == Res.best_param(curclass));
    
    Res.bestsigma(curclass) = imrelief.sigma(xpos{curclass}(1));
    Res.bestlambda(curclass) = imrelief.lambda(ypos{curclass}(1));
    
    fprintf('\nOptimum %g at sigma = %g, lambda = %g.', mxMxClx, ...
        Res.bestsigma(curclass), ...
        Res.bestlambda(curclass))
    
    if FEATSEL.salthreshmode
        Res.bestthresh(curclass) = Res.mThreshCl{curclass}(xpos{curclass}(1),ypos{curclass}(1),bestd);
        Res.bestcrit(curclass) = Res.mCritCl{curclass}(xpos{curclass}(1),ypos{curclass}(1),bestd);
        fprintf('\nAverage threshold at optimum = %g (Crit = %g)',Res.bestthresh(curclass),Res.bestcrit(curclass))
    end
    
    % Display grids
    
    StrIntro = ['Smoothed IMRelief grid of classifier [ ' CV.class{1,1}{curclass}.groupdesc ...
        ' ] \newlineMax = ' num2str(Res.best_param(curclass),'%1.2f') ...
        ' at [ lambda = ' num2str(Res.bestlambda(curclass)) ...
        ', sigma = ' num2str(Res.bestsigma(curclass)) ];
    
    if ndx > 1
        StrIntro = [StrIntro ', dim = ' num2str(dx(bestd)) ' ]'];
    else
        StrIntro = [StrIntro ' ]'];
    end
    
    rangez = [min(Res.best_MxClx{curclass}(:)) max(Res.best_MxClx{curclass}(:))];
    try
        display_grid(h, nclass, curclass, StrIntro, rangez, ...
            imrelief.sigma, ...
            imrelief.lambda, ...
            Res.best_MxClx{curclass})
    catch
       fprintf('\nPrinting failed')
    end
    
end

end
% _________________________________________________________________________

function display_grid(fignum, nclass, curclass, strintro, rangez, Cs, Gs, grid)

set(0,'CurrentFigure',fignum);
%subplot(1,nclass,curclass), contourf(1:numel(Gs),1:numel(Cs),grid)
subplot(1,nclass,curclass), bar3color(grid,1)
set(gca,'XTickLabel',Gs);
%set(gca,'XDir','reverse','XTick',Gs)
set(gca,'YTickLabel',Cs);
%set(gca,'YLim',[Cs(1) Cs(end)])
pos = get(gca,'Position');
pos(4)=pos(4)-0.08;
set(gca,'Position',pos);
ti=get(gca,'TightInset');
set(gca,'LooseInset',ti+0.02)
title(strintro, 'HorizontalAlignment','center','FontSize',10);
xlabel('lambda','FontSize',10);
ylabel('sigma','FontSize',10);
if rangez(2) > rangez(1)
    set(gca,'CLim',[rangez(1),rangez(2)])
end
xlim([.5 numel(Gs)+.5])
ylim([.5 numel(Cs)+.5])
%clim([rangez(1) rangez(2)])
colorbar('Location','SouthOutside');
colormap('autumn')

end

function C = compute_weightvecsim(D1,D2)

a = size(D1); a(1)=[];
C = zeros(a(1),a(2));
for i=1:prod(a)
    C(i) = pdist([D1(:,i)';D2(:,i)'],'cosine');
end

end