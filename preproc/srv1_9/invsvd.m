function Kinv=invsvd(K)
% compute matrix inverse of K using SVD
% Yifeng Li

[U,S,V]=svd(K);
Sinv=S;
[m,n]=size(S);
rk=m;
for i=1:min(m,n)
    if S(i,i)>eps
        Sinv(i,i)=1/S(i,i);
    else
        rk=i;
        break;
%         Sinv(i,i)=0;
    end
end
% Kinv=U*Sinv*V';
% rk=rank(K);
Kinv=U(:,1:rk)*Sinv(1:rk,1:rk)*V(:,1:rk)';
end