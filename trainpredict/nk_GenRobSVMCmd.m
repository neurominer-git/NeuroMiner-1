function cmd = nk_GenRobSVMCmd(ROBSVM, Params)

options = nk_GenMatLearnOptions(Params);

cmd = sprintf('%s %s', options.learner, options.kernelFunc);
switch options.learner
    case '-s 0'
        cmd = sprintf('-c %g %s', options.cost, cmd);
    case '-s 1'
        cmd = sprintf('-nu %g %s', options.nu, cmd);
    case '-s 3'
        cmd = sprintf('-c %g -p %g, %s', options.cost, options.epsilon, cmd);
    case '-s 4'
        cmd = sprintf('-c %g -nu %g, %s', options.cost, options.nu, cmd);
end
switch options.kernelFunc
    case '-t 1'
        cmd = sprintf('%s -r %g -d %g', options.coef0, options.degree);
    case '-t 2'
        cmd = sprintf('%s -g %g', options.gamma);
end

cmd = sprintf('%s -m %g -e %g -h %g -b %g', cmd, ...
                        ROBSVM.Optimization.m, ...
                        ROBSVM.Optimization.e, ...
                        ROBSVM.Optimization.h, ...
                        ROBSVM.Optimization.b);