    % Non-negative Matrix Factorization via Nesterov's Optimal Gradient Method.
    % NeNMF: Matlab Code for Efficient NMF Solver

    % Reference
    %  N. Guan, D. Tao, Z. Luo, and B. Yuan, "NeNMF: An Optimal Gradient Method
    %  for Non-negative Matrix Factorization", IEEE Transactions on Signal
    %  Processing, Vol. 60, No. 6, PP. 2882-2898, Jun. 2012. (DOI:
    %  10.1109/TSP.2012.2190406)

    %   Modified by:F. Yahaya
    %   Date: 06/09/2018
    %   Contact: farouk.yahaya@univ-littoral.fr


    % <Inputs>
    %        X : Input data matrix (m x n)
    %        r : Target low-rank
    %
    %        MAX_ITER : Maximum number of iterations. Default is 1,000.
    %        MIN_ITER : Minimum number of iterations. Default is 10.

    %        TOL : Stopping tolerance. Default is 1e-5. If you want to obtain a more accurate solution, decrease TOL and increase MAX_ITER at the same time.

    % <Outputs>
    %        W : Obtained basis matrix (m x r).
    %        H : Obtained coefficients matrix (r x n).
    %        T : CPU TIME.
    %        RRE: Relative reconstruction error in each iteration


    %        Tmax : CPU time in seconds.
    % Note: another file 'stop_rule.m' should be included under the same
    % directory as this code.




    function [W,H,RRE,T]=VANILLA_NeNMF(X, W, H, r, Tmax)
    MinIter=10;

    tol=1e-5;
    T=zeros(1,301);
    RRE=zeros(1,301);
    
    if isempty(W) || isempty(H)
        [m,n] = size(X);
        rng(1)
        W = rand(m,r);
        H = rand(r,n);
    end

    ITER_MAX=500;      % maximum inner iteration number (Default)
    ITER_MIN=10;        % minimum inner iteration number (Default)

    HVt=H*X'; HHt=H*H';
    WtV=W'*X; WtW=W'*W;

    GradW=W*HHt-HVt';
    GradH=WtW*H-WtV;

    init_delta=stop_rule([W',H],[GradW',GradH]);
    tolH=max(tol,1e-3)*init_delta;
    tolW=tolH;% Stopping tolerance


    W=W';
    k=1;
    tic
    RRE(k) = nmf_norm_fro( X, W', H);
    T(k) = 0;
    % main loop
    while(toc<= Tmax+0.05)
          
        % Optimize H with W fixed
        [H,iterH]=NNLS(H,WtW,WtV,ITER_MIN,ITER_MAX,tolH);
        
        if iterH<=ITER_MIN
            tolH=tolH/10;
        end
        
        HHt=H*H';   HVt=H*X';
        
        % Optimize W with H fixed
        [W,iterW,GradW]=NNLS(W,HHt,HVt,ITER_MIN,ITER_MAX,tolW);
        if iterW<=ITER_MIN
            tolW=tolW/10;
        end
        WtW=W*W'; WtV=W*X;
        GradH=WtW*H-WtV;
        %     HIS.niter=niter+iterH+iterW;
        delta=stop_rule([W,H],[GradW,GradH]);
        % Output running detials
        
        %     % Stopping condition
        if (delta<=tol*init_delta && k>=MinIter)
            break;
        end

        if toc-(k-1)*0.05>=0.05
            k = k+1;
            RRE(k) = nmf_norm_fro( X, W', H);
            T(k) = toc;
        end
        
    end  %end of  loop


    W=W';
    return;

    function [H,iter,Grad]=NNLS(Z,WtW,WtV,iterMin,iterMax,tol)


    if ~issparse(WtW)
        L=norm(WtW);	% Lipschitz constant
    else
        L=norm(full(WtW));
    end
    H=Z;    % Initialization
    Grad=WtW*Z-WtV;     % Gradient
    alpha1=1;

    for iter=1:iterMax
        H0=H;
        H=max(Z-Grad/L,0);    % Calculate sequence 'Y'
        alpha2=0.5*(1+sqrt(1+4*alpha1^2));
        Z=H+((alpha1-1)/alpha2)*(H-H0);
        alpha1=alpha2;
        Grad=WtW*Z-WtV;
        
        % Stopping criteria
        if iter>=iterMin
            % Lin's stopping condition
            pgn=stop_rule(Z,Grad);
            if pgn<=tol
                break;
            end
        end
    end

    Grad=WtW*H-WtV;

    return;
    function f = nmf_norm_fro(X, W, H)
    % Author : F. Yahaya
    % Date: 13/04/2018
    % Contact: farouk.yahaya@univ-littoral.fr
    % Goal: compute a normalized error reconstruction of the mixing matrix V
    % "Normalized" means that we divide the squared Frobenius norm of the error
    % by the squared Frobenius norm of the matrix V
    % Note: To express the error in dB, you have to compute 10*log10(f)
    %

    f = norm(X - W * H,'fro')^2/norm(X,'fro')^2;

    return;
