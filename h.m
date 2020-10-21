function hX = h(rX, epsilon)

lenVec=size(rX,1);
numVec=size(rX,2);
% Compute the conjunction of all binary
% rank distance matrices.
Cx=0;
for n=1:lenVec
    xBinVec = abs(rX(1:lenVec-n,1)-rX(1+n:lenVec,1))<epsilon;
    for i = 2:numVec
        xBinVec = xBinVec & (abs(rX(1:lenVec-n,i)-rX(1+n:lenVec,i))<epsilon);
    end 
    % Successively update correlation integral. 
    Cx=Cx+sum(xBinVec);
end 

% Finally, compute the entropy from the value 
% of the correlation integral Cx.

hX = -log2((1.0/lenVec)*(1.0+((2.0/lenVec)*Cx)));

return