function [G H w] = mnmf(X,r,labl,C,K)

theta = C/K;
str = ['-t', ' ', '0', ' ', '-c', ' ', num2str(C)];
fprintf('K = %f     C = %f     theta = %f\n',K,C,theta);
MAX_ITER = 30;
[siz_train_dim siz_no_fea] = size(X);
Ginit = rand(siz_train_dim,r);
Hinit = rand(r,siz_no_fea);
[m,n] = size(X);
G = Ginit;
H = Hinit;

% [aa bb] = size(H*H');
for iter = 1:MAX_ITER

%<------------------- Solving for G -------------------
    if(iter==1)
%         G = X*H'*inv(H*H'+0.0000000001*eye(aa));
        G = X*H'*inv(H*H');
    else
%         B        =       inv(2*H*H'+0.0000000001*eye(aa));
        B        =       inv(2*H*H');
        E        =       2*X*H';
        HB       =       H'*B;
        HBw      =       HB*w;
        XX       =       X'*X;
        BE       =       B*E';
        BEX      =       BE*X;
        HBEX     =       H'*BEX;
        BH       =       HB';
        wBw      =       w'*B*w;
        wB       =       w'*B;
        wBE      =       wB*E';
        wBH      =       w'*BH;
        sq_wBH   =       wBH*wBH';
        

        vec1 = zeros(n,1);
        vec2 = vec1;
        vec3 = vec1;
        
        
        for i = 1:n
            vec1(i,1) = labl(1,i)*sum(XX(i,:)'.*HBw);
            vec2(i,1) = labl(1,i)*sum(HBEX(:,i).*wBH');
            vec3(i,1) = labl(1,i)*wBE*X(:,i);
        end
        vec4 = b*labl';
        vec5 = ones(n,1);
        vec1 = 2*vec1;
        vec2 = 2*vec2;

        
        Hq_temp = zeros(n);
        for i = 1:n
            for j = 1:n
                Hq_temp(i,j) = -(labl(1,i)*labl(1,j)*XX(i,j)*sq_wBH) + (labl(1,i)*labl(1,j)*...
                    wBw*XX(i,j));
            end
        end
        fq = vec1 - vec2 + vec3 + vec4 - vec5;
        Hq = 2*Hq_temp;
        clear Hq_temp;
        lq = zeros(n,1);
        uq = theta*ones(n,1);
        [alpha,err,lm] = qpas(Hq,fq,[],[],[],[],lq,uq,0);
        alpha = max(alpha,eps);
        temp_F = zeros(m,r);
        for i = 1:n
            temp_F = temp_F + labl(1,i)*alpha(i,1)*X(:,i)*w';
        end
        G = (2*X*H' + temp_F)*B;
    end
%------------------- Solving for G ------------------->


%<------------------ Solving for w,b ---------------
    % w: svm hyperplane vector, r x 1
    % b: svm parameter, 1 x 1
    mat_svm = G'*X;
    [w b] = calling_svm(mat_svm',labl',str);
%----------------- Solving for w,b,psi ----------------->

%<------------------ Solving for H ---------------------
    B = inv(2*G'*G);
    Hq = (B+B')/2;
    Fq = 2*B*G'*X;
    temp1 = 2*G'*X;
    for i = 1:n
        fq = Fq(:,i);
        lq = zeros(r,1);
        [gamma,err,lm] = qpas(Hq,fq,[],[],[],[],lq,[],0);
        gamma = max(gamma,eps);
        H(:,i) = B*(temp1(:,i) + gamma);
    end
    H = max(H,eps); 
%------------------ Solving for H --------------------->   
fprintf('.');
end
fprintf('\n');



%------------------ function calling_svm() --------------------->   

function [w b] = calling_svm(H,labl,str)
model = svmtrain291(labl,H,str);    
w = (model.SVs' * model.sv_coef);
b = -model.rho;
if model.Label(1) == -1
  w = -w;
  b = -b;
end
