function [Y,mindat,maxdat] = scaledata(Y,cv,minval,maxval,mindat,maxdat)
%
% Program to scale the values of a matrix from a user specified minimum to a user specified maximum.
% Optional scaling to the mix-max range of a different matrix using the parameter mindat / maxdat of that matrix
%
% Usage:
% [outputData,[mindat],[maxdat]] = scaleData(inputData,minVal,maxVal, [mindat], [maxdat]);
%
% Example:
% a = [1 2 3 4 5];
% a_out = scaledata(a,0,1);
% 
% Output obtained: 
%            0    0.1111    0.2222    0.3333    0.4444
%       0.5556    0.6667    0.7778    0.8889    1.0000
%
% Program written by:
% Aniruddha Kembhavi, July 11, 2007,
% Modified by Nikos Koutsouleris, May 23, 2008

if isstruct(Y)
    if nargin < 6,maxdat = max(Y.mapY(:));end;
    if nargin < 5,mindat = min(Y.mapY(:));end;
    [Y.mapY,mindat,maxdat] = scaleY(Y.mapY,minval,maxval,mindat,maxdat);
    nperms=size(Y.Y,1);
    fold=size(Y.Y,2);
    for i=1:nperms
        for j=1:fold
            [Y.Y{i,j}(cv.TrainInd{i,j},:),mndat,mxdat] = scaleY(Y.Y{i,j}(cv.TrainInd{i,j},:),minval,maxval,mindat,maxdat);
            Y.Y{i,j}(cv.TestInd{i,j},:) = scaleY(Y.Y{i,j}(cv.TestInd{i,j},:),minval,maxval,mndat,mxdat);
        end
    end
else
    if nargin < 6,maxdat = max(Y(:));end;
    if nargin < 5,mindat = min(Y(:));end;
    [Y,mindat,maxdat] = scaleY(Y,minval,maxval,mindat,maxdat);
end

return

function [dataout,mindat,maxdat] = scaleY(datain,minval,maxval,mindat,maxdat)

rangedat = maxdat- mindat;
scalefactor = repmat(maxval-minval,1,size(rangedat,2))./rangedat;

dataout = datain - mindat;
dataout = dataout * scalefactor;
dataout = dataout + minval;

return