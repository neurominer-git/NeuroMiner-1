function p_map = nk_ComputeSVMPmap( Y, labels, model)

p       = sum(labels==1)/max(size(labels));
if ~isfield(model,'w')
    w   = model.SVs' * model.sv_coef;
else
    w   = model.w';
end
X       = Y;
[r,c]   = size(Y);
J       = ones(r,1);
K       = X*X';
Z       = inv(K)+(inv(K)*J*inv(-J'*inv(K)*J)*J'*inv(K));
C       = X'*Z;
SD      = sqrt(sum(C.^2,2)*(4*p-4*p^2));
mean_w  = sum(C,2)*(2*p-1);
p_map   = 2*normcdf(-abs(w-mean_w),zeros(size(w)),SD);

end