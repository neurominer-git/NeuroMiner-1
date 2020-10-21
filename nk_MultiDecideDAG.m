function [P, Perf, Sim] = nk_MultiDecideDAG(H, Y, Classes, ngroups)
% function [P, Perf, d, crit] = nk_MultiDecideOneVsOne(H, Y, Classes,
% maxfunc, weightflag)
%
% Perform One-Vs-One multi-group classification
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Nikolaos Koutsouleris, 01/2019
global SVM

if iscell(H)
    [ iy , jy ]     = size(H);
    P               = cell(iy,jy);
    Sim             = cell(iy,jy);
    Perf            = zeros(iy,jy);
else
    iy = 1; jy=1;
end

if ~isempty(SVM) && SVM.GridParam == 14
    multimode = 1;
else
    multimode = 0;
end

for k=1:iy % Loop through partitions (in case iscell(H) = true)
    
    for l=1:jy % Loop through folds (in case iscell(H) = true)
        
        if iscell(H)
            Hkl = H{k,l};
            Ckl = Classes{k,l};
        else
            Hkl = H;
            Ckl = Classes;
        end
        
        % Construct coding matrix
        [~, T] = nk_OneVsOne(Ckl,ngroups);
        nsubj       = size(Hkl,1);
                
        % generate graph
        G = generateDAG(T); nG = size(G,2);
        td = ones(nsubj,nG);
        
        for i = 1 : nG-1 
            
            di = G(:,i);    divec = find(di);
            dj = G(:,i+1);  djvec = find(dj);
            cnt = 1;
            for j=1:numel(divec)
                I = Ckl == di(divec(j));
                Cj = sign(sum(Hkl(:,I),2)); 
                indp = Cj > 0 & td(:,i) ==di(divec(j)); 
                indn = Cj < 0 & td(:,i) ==di(divec(j)); 
                td(indp,i+1) = dj(djvec(cnt)); cnt=cnt+1;
                td(indn,i+1) = dj(djvec(cnt)); cnt=cnt+1;
            end
        end
        
        tdx = td(:,nG); utdx = unique(tdx); nutdx = numel(utdx);
        Px = zeros(size(tdx)); SimX = zeros(size(tdx,1), ngroups);
        for i=1:nutdx
            
            % Get current dichotomizer and classes it decides upon
            Mi = T(:,utdx(i));
            Cp = find(Mi==1);
            Cn = find(Mi==-1);
            
            % Index to applicable dichotomizers in coding matrix
            Cj = Ckl == utdx(i);
            
            % Index to applicable dichotomizers in sample
            I = tdx == utdx(i); 
            
            % Voting 
            Cji = sign(sum(Hkl(:,Cj),2));
            
            % Assign classes
            Px(I & Cji==1) = Cp; 
            Px(I & Cji==-1) = Cn;
            SimX(I & Cji == 1, Cp) = 1;
            SimX(I & Cji ==-1, Cn) = 1;
        end
        
        if iscell(H)
            P{k,l} = Px;
            Sim{k,l} = SimX;
            if isempty(Y), Y = ones(numel(P),1); end
            Perf(k,l) = nk_MultiPerfQuant(Y, P, multimode);
        else
            P = Px;
            Sim = SimX;
            if isempty(Y), Y = ones(numel(P),1); end
            Perf = nk_MultiPerfQuant(Y, P, multimode);
        end
    end
end

