function GridAct = nk_GenGridAct_batch(cv, curCPU, numCPU, x1, x2, y1, y2, xx1, xx2, yy1, yy2)

% %%%%%% DEFINE CV ACTION GRID ACCORDING TO "numCPU" and "curCPU" %%%%%%%%

if nargin < 8 
    % Outer CV mode
    [ix, jx] = size(cv.TrainInd);
    ixP = x2-(x1-1); jxP = y2-(y1-1);
    numCV2parts = ixP * jxP;
    numCV2partsCPU = ceil(numCV2parts / numCPU);
    GridActP = false(numCV2parts,1);
    startpos = 1 + (curCPU-1) * numCV2partsCPU;
    endpos = curCPU * numCV2partsCPU;
    if endpos > numCV2parts, endpos = numCV2parts; end
    GridActP( startpos : endpos) = true;
    GridActP = reshape(GridActP,ixP,jxP);
    GridAct = false(ix,jx);
    GridAct(x1:x2,y1:y2) = GridActP;
else % Inner CV mode
    [ix, jx] = size(cv.TrainInd);
    [iy, jy] = size(cv.cvin{1,1}.TrainInd);
    ixP = x2-(x1-1); jxP = y2-(y1-1);
    iyP = xx2-(xx1-1); jyP = yy2-(yy1-1);
    numCVparts = ixP * jxP * iyP * jyP;
    numCVpartsCPU = ceil(numCVparts / numCPU);
    GridActP = false(numCVparts,1);
    startpos = 1 + (curCPU-1) * numCVpartsCPU;
    endpos = curCPU * numCVpartsCPU;
    if endpos > numCVparts, endpos = numCVparts; end
    GridActP( startpos : endpos) = true;
    GridActP = reshape(GridActP, ixP, jxP, iyP, jyP);
    GridAct = false(ix, jx, iy, jy);
    GridAct(x1:x2, y1:y2, xx1:xx2, yy1:yy2) = GridActP;
end
return
