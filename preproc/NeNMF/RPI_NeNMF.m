    % Randomized Power Iterations NeNMF (RPINeNMF)

    % Reference
    %  N. Guan, D. Tao, Z. Luo, and B. Yuan, "NeNMF: An Optimal Gradient Method
    %  for Non-negative Matrix Factorization", IEEE Transactions on Signal
    %  Processing, Vol. 60, No. 6, PP. 2882-2898, Jun. 2012. (DOI:
    %  10.1109/TSP.2012.2190406)

    %   Modified by:F. Yahaya
    %   Date: 06/09/2018
    %   Contact: farouk.yahaya@univ-littoral.fr


    %  Reference
    %  F. Yahaya, M. Puigt, G. Delmaire, G. Roussel, Faster-than-fast NMF using
    %  random projections and Nesterov iterations, to appear in the Proceedings
    %  of iTWIST: international Traveling Workshop on Interactions between
    %  low-complexity data models and Sensing Techniques, Marseille, France,
    %  November 21-23, 2018


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

    function [W,H,RRE,T]=RPI_NeNMF( X,W,H,r, Tmax)
    MinIter=10;

    % tol=1e-5;
    tol = 1e-5;
    if isempty(W) || isempty(H)
        [m,n] = size(X);
        rng(200)
        W = rand(m,r);
        H = rand(r,n);
    end
    T=zeros(1,5000);
    RRE=zeros(1,5000);


    ITER_MAX=500;      % maximum inner iteration number (Default)
    ITER_MIN=10;        % minimum inner iteration number (Default)

    [L,R]=compression(X,r);

    X_L = L * X;
    X_R = X* R;


    %initialization
    % W=W0; H=H0;

    H_comp= H* R;
    W_comp = L*W;

    HVt=H_comp*X_R';
    HHt=H_comp*H_comp';

    WtV=W_comp'*X_L;
    WtW=W_comp'*W_comp;

    GradW=W*HHt-HVt';
    GradH=WtW*H-WtV;

    init_delta=stop_rule([W',H],[GradW',GradH]);
    tolH=max(tol,1e-3)*init_delta;
    tolW=tolH;               % Stopping tolerance


    % Iterative updating


    W=W';
    k=1;
    RRE(k) = nmf_norm_fro( X, W', H);
    T(k) =0;
    tic
    % main loop
    while (toc<= Tmax)
        
        k = k+1;
        
        % Optimize H with W fixed
        [H,iterH]=NNLS(H,WtW,WtV,ITER_MIN,ITER_MAX,tolH);
        
        if iterH<=ITER_MIN
            tolH=tolH/10;
        end
        H_comp=H*R;
        HHt=H_comp*H_comp';   HVt=H_comp*X_R';
        % Optimize W with H fixed
        [W,iterW,GradW]=NNLS(W,HHt,HVt,ITER_MIN,ITER_MAX,tolW);
        
        if iterW<=ITER_MIN
            tolW=tolW/10;
        end
        W_comp=W * L';
        WtW=W_comp*W_comp';
        WtV=W_comp*X_L;
        GradH=WtW*H-WtV;
        
        %     HIS.niter=niter+iterH+iterW;
        delta=stop_rule([W,H],[GradW,GradH]);
        
        % Stopping condition
        if (delta<=tol*init_delta && k>=MinIter)
            break;
        end
        
        
        RRE(k) = nmf_norm_fro( X, W', H);
        T(k) = toc;
        
        
        
    end  %end of  loop
    W=W';


    end


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

    end


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

    end

    function [ L,R ] = compression(X, r)
    %             Tepper, M., & Sapiro, G. (2016). Compressed nonnegative
    %             matrix factorization is fast and accurate. IEEE Transactions
    %             on Signal Processing, 64(9), 2269-2283.
    %             see: https://arxiv.org/pdf/1505.04650
    %             The corresponding code is originally created by the authors
    %             Then, it is modified by F. Yahaya.
    %             Date: 13/04/2018
    %

    compressionLevel=20;
    [m,n]=size(X);

    l = min(n, max(compressionLevel, r + 10));

    OmegaL = randn(n,l);

    B = X * OmegaL;
    for j = 1:4
        B = X * (X' * B);
    end
    [L,~] = qr(B, 0);
    L = L';

    OmegaR = randn(l, m);
    B = OmegaR * X;
    for j = 1:4
        B = (B * X') * X;
    end

    [R,~] = qr(B', 0);
    end
