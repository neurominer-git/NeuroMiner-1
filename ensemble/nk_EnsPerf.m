function [hEPerf, hE] = nk_EnsPerf(E, L)
global MODEFL EVALFUNC

% Compute ensemble performance
switch MODEFL
    case 'classification'
        hE = sign(mean(E,2));
        % throw a coin on the zeros
        hE = nk_ThrowCoin(hE);

    case 'regression'
        hE = mean(E,2);
       
end

% Measure accuracy or some other criterion
hEPerf = feval(EVALFUNC, L, hE);

return