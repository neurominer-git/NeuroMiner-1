% Given new scores yHat, and a mapping yHatIR -> yIR, calibrates the scores
% in yHat
function yCal = calibrateScores(yHat, yHatIR, yIR)

    n = length(yHatIR);

    [v,i] = sort(yHatIR);
    yHatIR = yHatIR(i); yIR = yIR(i);
    
    yCalIdx = bsearch(yHatIR, yHat);

    delta = yHat - yHatIR(yCalIdx); %disp('delta:'); disp(delta(1:10));
    leftSide = delta < 0;
    rightSide = delta >= 0;

    % Scale the delta's to represent the fraction within in the appropriate
    % interval
    delta(leftSide) = abs(delta(leftSide)) ./ (yHatIR(yCalIdx(leftSide)) - yHatIR(max(1,yCalIdx(leftSide) - 1)));
    delta(rightSide) = abs(delta(rightSide)) ./ (yHatIR(min(n, yCalIdx(rightSide) + 1)) - yHatIR(yCalIdx(rightSide)));

    delta(isnan(delta)) = 0;
    delta(isinf(delta)) = 0;

    yCal = zeros(size(yHat));
    yCal(leftSide) = delta(leftSide) .* yIR(yCalIdx(leftSide)) + (1 - delta(leftSide)) .* yIR(max(1,yCalIdx(leftSide) - 1));   
    yCal(rightSide) = delta(rightSide) .* yIR(yCalIdx(rightSide)) + (1 - delta(rightSide)) .* yIR(min(n, yCalIdx(rightSide) + 1));    
