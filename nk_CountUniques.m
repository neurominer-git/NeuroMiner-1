function R = nk_CountUniques(Y, PercThresh)

N = isfinite(Y);
U = zeros(size(Y,2),1); UX = cell(size(Y,2),1); UP = cell(size(Y,2),1);
for i=1:size(Y,2), 
    Yi = Y(N(:,i),i);
    UX{i} = unique(Yi);
    UN{i}= zeros(size(UX{i}));
    UP{i}= zeros(size(UX{i}));
    for j=1:numel(UX{i})
        UN{i}(j) =  sum(Yi==UX{i}(j));
        UP{i}(j) =  UN{i}(j)*100/numel(Yi);
        %UP{i}(j) = sum(Yi==UX{i}(j));
    end
    U(i)  = numel(UX{i});
end

R.U = U;
R.UX = UX;
R.UN = UN;
R.UP = UP;

if exist('PercThresh','var') && ~isempty(PercThresh)
    UT = false(size(Y,2),1);
    for i=1:size(Y,2)
        if sum(UP{i} >= PercThresh) 
            UT(i) = true; 
        end
    end
    R.UT = UT;
end
    