function param = nk_PolyKernel_config(param)

switch param.kernel.kernstr 
    case {' -t 1', 'poly', 'polynomial', 'Polynomial', 'hpolyN', 'polyN'}
        param.kernel.poly_coef = nk_input('Define coefficient for polynomial kernel',0,'e',0);
end

end