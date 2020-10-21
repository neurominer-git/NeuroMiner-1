function [ sY, IN ] = PerfSymbolizeObj(Y, IN)

if isempty(IN),eIN=true; else, eIN=false; end
% Default params for binning
if eIN|| ~isfield(IN,'SYMBOL') || isempty(IN.SYMBOL),  
    IN.SYMBOL.symMinBin = 3;
    IN.SYMBOL.symMaxBin  = 25;
    IN.SYMBOL.symSeqLength  = 4;
    IN.SYMBOL.symStdNum  = 3;
end

if eIN || ~isfield(IN,'sBin') || isempty(IN.sBin)
    [ sY, IN.sBin ] = symbolize(Y, IN.SYMBOL.symMinBin, IN.SYMBOL.symMaxBin, IN.SYMBOL.symSeqLength, IN.SYMBOL.symStdNum);
else
    sY = zeros(size(Y));
    for j=1:size(Y,2)
        binvec = IN.sBin(1,j) : IN.sBin(2,j) : IN.sBin(3,j);
        binvecx = numel(binvec)-1;
        if numel(binvec)~=0
            % Assign values below binvec(1) to 1
			sY(Y(:,j) < binvec(1),j) = 1;
			for k=1 : binvecx, sY(Y(:,j) >= binvec(k) & Y(:,j) < binvec(k+1),j ) = k; end
			% Assign values above binvec(k) to binvec(k)
			sY(Y(:,j) >= binvec(k+1),j) = k; 
        end
    end
end
