function [eigvector, elapse] = KSR(options,Responses,K)
% KSR: Kernel Spectral Regression
%
%       [eigvector, elapse] = KSR(options,Responses,K)
% 
%             Input:
%                  K    - Kernel matrix. 
%           Responses   - response vectors. Each column is a response vector 
%
%               options - Struct value in Matlab. The fields in options
%                         that can be set:
%
%                          ReguAlpha   -  The regularization parameter.
%                                           Default value is 0.1.  
%
%                            ReguType  -  'Ridge': Tikhonov regularization
%                                                  L2-norm regularizer
%                                          Only support 'Ridge' now.
%
%
%
%             Output:
%               eigvector - Each column is an embedding function, for a new
%                           data point (row vector) x,  y = K(x,:)*eigvector
%                           will be the embedding result of x.
%                           K(x,:) = [K(x1,x),K(x2,x),...K(xm,x)]
%               elapse    - Time spent on different steps 
%                           
% 
%
%    Examples:
%
%
%
%Reference:
%
%
%   Deng Cai, Xiaofei He, Jiawei Han, "Efficient Kernel Discriminant
%   Analysis via Spectral Regression", Department of Computer Science
%   Technical Report No. 2888, University of Illinois at Urbana-Champaign
%   (UIUCDCS-R-2007-2888), August 2007.  
% 
%   Deng Cai, Xiaofei He, Jiawei Han, "Spectral Regression for
%   Dimensionality Reduction", Department of Computer Science
%   Technical Report No. 2856, University of Illinois at Urbana-Champaign
%   (UIUCDCS-R-2007-2856), May 2007.  
%
%   Deng Cai, Xiaofei He, Jiawei Han, "Semi-Supervised Regression using
%   Spectral Techniques", Department of Computer Science
%   Technical Report No. 2749, University of Illinois at Urbana-Champaign
%   (UIUCDCS-R-2007-2749), July 2006.  
%
%
%   version 2.0 --Aug/2007 
%   version 1.0 --May/2006 
%
%   Written by Deng Cai (dengcai2 AT cs.uiuc.edu)
%

bGlobalK = 0;
if ~exist('K','var')
    global K;
    bGlobalK = 1;
end

if ~exist('Responses','var')
    global Responses;
end


if ~isfield(options,'ReguAlpha')
    options.ReguAlpha = 0.001;
end

LassoCardi = 0;
if isfield(options,'LassoCardi') 
    LassoCardi = max(options.LassoCardi)+10;
end

if strcmpi(options.ReguType,'Lasso') 
    options.ReguAlpha = 0;
    options.ReguType = 'RidgeLasso';
end



elapse.timeChol = 0;
elapse.timeBack = 0;

tmp_T = cputime;

nVector = size(Responses,2);

switch lower(options.ReguType)
    case {lower('Ridge')}

        if options.ReguAlpha > 0
            for i=1:size(K,1)
                K(i,i) = K(i,i) + options.ReguAlpha;
            end
        end

        R = chol(K);
        elapse.timeChol = cputime - tmp_T;

        tmp_T = cputime;
        eigvector = R\(R'\Responses);
        elapse.timeBack = cputime - tmp_T;

        tmp_T = cputime;
        tmpNorm = sqrt(sum((eigvector'*K).*eigvector',2));
        eigvector = eigvector./repmat(tmpNorm',size(eigvector,1),1);
    case {lower('RidgeLasso')}
        if isfield(options, 'ResponseStart')
            iStart = options.ResponseStart;
        else
            iStart = 1;
        end
        if isfield(options, 'ResponseEnd')
            iEnd = options.ResponseEnd;
        else
            iEnd = nVector;
        end
        if isfield(options, 'ResponseRange')
            ResponseRange = options.ResponseRange;
        else
            ResponseRange = iStart:iEnd;
        end

        K = double(K);
        if options.ReguAlpha > 0
            for i=1:size(K,1)
                K(i,i) = K(i,i) + options.ReguAlpha;
            end
        end

        global Gram;
        Gram = K'*K;
        Gram = max(Gram,Gram');
%         K = single(K);

        eigvector = {};
        global X;
        X = K;
        if bGlobalK
            clear global K;
        else
            clear K;
        end
        
        for i = ResponseRange
            disp(['Total projections: ',num2str(nVector),'; Projection now: ',num2str(i)]);
            eigvector_T = lars([], Responses(:,i),'lasso', -LassoCardi,1,[],options.LassoCardi);
            if isfield(options,'tmpSave') & options.tmpSave
                save(['eigvector_',num2str(i),'.mat'],'eigvector_T');
            else
                eigvector{i} = eigvector_T;
            end
            clear eigvector_T;
        end
        clear global Gram;
        
        if bGlobalK
            global K;
            K = X;
        else
            K = X;
        end
        clear global X;

        if isfield(options, 'tmpExit') & options.tmpExit
            error('tmpExit');
        end

        if isfield(options,'tmpSave') & options.tmpSave
            for i = 1:nVector
                load(['eigvector_',num2str(i),'.mat'],'eigvector_T');
                eigvector{i} = eigvector_T;
            end
        end
    otherwise
        error('ReguType does not exist!');
end
        

if strcmpi(options.ReguType,'Lasso') | strcmpi(options.ReguType,'RidgeLasso')
    eigvectorAll = eigvector;
    eigvector = cell(length(options.LassoCardi),1);
    
    for i = 1:length(eigvectorAll)
        eigvector_T = full(eigvectorAll{i});
        [tm,tn] = size(eigvector_T);
        tCar = zeros(tn,1);
        for k = 1:tn
            tCar(k) = length(find(eigvector_T(:,k)));
        end

        for cardidx = 1:length(options.LassoCardi)
            ratio = options.LassoCardi(cardidx);
            iMin = find(tCar == ratio);
            if length(iMin) == 0
                error('Card dose not exist!');
            end
            tmpEigvec = eigvector_T(:,iMin(end))/sqrt(eigvector_T(:,iMin(end))'*K*eigvector_T(:,iMin(end)));
            eigvector{cardidx} = [eigvector{cardidx} tmpEigvec];
        end
    end
end


elapse.timeAll =  cputime - tmp_T + elapse.timeChol +  elapse.timeBack ;




