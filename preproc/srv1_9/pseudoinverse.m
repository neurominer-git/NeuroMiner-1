function Hinv=pseudoinverse(H,lambda)
% compute pseudoinverse of matrix H
% Usage:
% Hinv=pseudoinverse(H);
% Hinv=pseudoinverse(H,lambda);
% H: matrix
% lambda: scalar, the parameter to reach stable result and avoid
% singularity, its value should be very large, i.e. 2^32.
% Hinv, the pseudoinverse of H
% Yifeng Li

[r,c]=size(H);
if nargin<2
    if r<c
        Hinv=(H'*H)\H';
    else
        Hinv=H'/(H'*H);
    end
    
    return;
end
if r<c
    Hinv=((1/lambda).*eye(c,c)+H'*H)\H';
else
    Hinv=H'/((1/lambda).*eye(r,r)+H'*H);
end
end