function [ sY, IN ] = nk_PerfZeroOut(sY, IN)
% Zeros out completely non-finite features
global VERBOSE

mY = size(sY,1);

% Zero-out non-finite data
if IN.zerooutflag == 1
   if ~isfield(IN,'indnonfin')
        IN.indnonfin = sum(~isfinite(sY)) == mY;
   end
   if sum(IN.indnonfin)>0
        if VERBOSE, fprintf(' zero-out completely non-finite features... '); end
        zerofeat = zeros(mY,sum(IN.indnonfin));
        sY(:,IN.indnonfin) = zerofeat;
   end
end