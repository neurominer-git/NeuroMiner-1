function [catarray, numelems] = nk_cellcat(arr, arr2, dim, F)

[ix,jx] = size(arr);
if ~exist('dim','var') || isempty(dim), dim = 1; end
catarray = [];
mx = 0; my = 0; cx = 0; cy = 0;
sx = zeros(ix*jx,1); sy= zeros(ix*jx,1);
% Determine output sizes of catarray
cnt=1;
for i=1:ix
    for j=1:jx
        if ~exist('arr2','var') || isempty(arr2)
            if ~exist('F','var') || isempty(F)
                [sx(cnt), sy(cnt)] = size(arr{i,j});
            else
                [sx(cnt), sy(cnt)] = size(arr{i,j}(:,F{i,j}));
            end
        else
            [sx1, sy1] = size(arr{i,j});
            [sx2, sy2] = size(arr2{i,j});
            sx(cnt) = sx1; sy(cnt) = sy1+sy2;
        end
        if sx(cnt) > mx, mx = sx(cnt); end
        if sy(cnt) > my, my = sy(cnt); end
        cx = cx + sx(cnt);
        cy = cy + sy(cnt);
        cnt=cnt+1;
    end
end

% Initialize with NaNs
switch dim
    case 1
        catarray = nan(cx,my);
    case 2
        catarray = nan(mx,cy);
end
cnt=1; startpos = 1; endpos = 0;
try
    for i=1:ix
        for j=1:jx
            
            if ~exist('arr2','var') || isempty(arr2)
            
                if ~exist('F','var') || isempty(F)
                    
                    switch dim
                        case 2
                            endpos = endpos + sy(cnt);
                            catarray(1:sx(cnt),startpos:endpos) = arr{i,j};
                            startpos = startpos + sy(cnt);
                        case 1
                            endpos = endpos + sx(cnt);
                            catarray(startpos:endpos,1:sy(cnt)) = arr{i,j};
                            startpos = startpos + sx(cnt);
                    end
                else
                    switch dim
                        case 2
                            endpos = endpos + sy(cnt);
                            catarray(1:sx(cnt),startpos:endpos) = arr{i,j}(:,F{i,j});
                            startpos = startpos + sy(cnt);
                        case 1
                            endpos = endpos + sx(cnt);
                            catarray(startpos:endpos,1:sy(cnt)) =  arr{i,j}(:,F{i,j});
                            startpos = startpos + sx(cnt);
                    end
                end
            else
                switch dim
                    case 2
                        endpos = endpos + sy(cnt);
                        catarray(1:sx(cnt),startpos:endpos) = [arr{i,j}; arr2{i,j}];
                        startpos = startpos + sy(cnt);
                    case 1
                        endpos = endpos + sx(cnt);
                        catarray(startpos:endpos,1:sy(cnt)) = [arr{i,j}; arr2{i,j}];
                        startpos = startpos + sx(cnt);
                end
            end
            cnt = cnt+1;
        end
    end
    switch dim
        case 1
            numelems = size(catarray,1);
        case 2
            numelems = size(catarray,2);
    end
catch
     catarray = []; numelems = 0;  
end
return