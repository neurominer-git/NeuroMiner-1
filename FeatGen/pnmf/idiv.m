function d = idiv(p0,q0)
eps = 1e-9;
p = p0(:) + eps;
q = q0(:) + eps;
d = sum(p.*log(p))-sum(p.*log(q))-sum(p)+sum(q);