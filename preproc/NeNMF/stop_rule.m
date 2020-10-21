    function retVal=stop_rule(X,gradX)
    % Stopping Criterions
    % Written by Naiyang (ny.guan@gmail.com)


    pGrad=gradX(gradX<0|X>0);
    retVal=norm(pGrad);

    end
