function [eigvector, elapse, itn, Y] = SR(options, Responses, data)
% SR: Spectral Regression
%
%       [eigvector, elapse, itn, Y] = SR(options, Responses, data)
% 
%             Input:
%               data    - data matrix. Each row vector of data is a
%                         sample vector. 
%           Responses   - response vectors. Each column is a response vector 
%
%               options - Struct value in Matlab. The fields in options
%                         that can be set:
%
%                        ReguAlpha    -  The regularization parameter.
%                                           Default value is 0.1.  
%
%                            ReguType  -  'Ridge': Tikhonov regularization
%                                                  L2-norm regularizer     
%                                         'Lasso': L1-norm regularizer
%                                    'RidgeLasso': Combine Ridge and Lasso
%                                         'Custom': User provided
%                                                   regularization matrix
%                                          Default: 'Ridge'
%
%                                          'Lasso' and 'RidgeLasso' will produce
%                                          sparse solution [See Ref 8]
%
%                        regularizerR  -   (nFea x nFea) regularization
%                                          matrix which should be provided
%                                          if ReguType is 'Custom'. nFea is
%                                          the feature number of data
%                                          matrix
%
%
%
%             Output:
%               eigvector - Each column is an embedding function, for a new
%                           sample vector (row vector) x,  y = x*eigvector
%                           will be the embedding result of x.
%
%                           If 'Lasso' or 'RidgeLasso' regularization is
%                           used, the output eigvector will be a cell,
%                           please see 'lars.m' for the result format.
%
%               elapse    - Time spent on different steps 
%                           
% 
%
%    Examples:
%
%               See SR_caller.m
%
%Reference:
%
%   1. Deng Cai, Xiaofei He, Jiawei Han, "Semi-Supervised Regression using
%   Spectral Techniques", Department of Computer Science
%   Technical Report No. 2749, University of Illinois at Urbana-Champaign
%   (UIUCDCS-R-2007-2749), July 2006.  
%
%   2. Deng Cai, Xiaofei He, Jiawei Han, "Spectral Regression for
%   Dimensionality Reduction", Department of Computer Science
%   Technical Report No. 2856, University of Illinois at Urbana-Champaign
%   (UIUCDCS-R-2007-2856), May 2007.  
%
%   3. Deng Cai, Xiaofei He, Jiawei Han, "SRDA: An Efficient Algorithm for
%   Large Scale Discriminant Analysis", Department of Computer Science
%   Technical Report No. 2857, University of Illinois at Urbana-Champaign
%   (UIUCDCS-R-2007-2857), May 2007.  
%
%   4. Deng Cai, Xiaofei He, and Jiawei Han. "Isometric Projection", Proc.
%   22nd Conference on Artifical Intelligence (AAAI'07), Vancouver, Canada,
%   July 2007.  
%
%   5. Deng Cai, Xiaofei He, Jiawei Han, "Efficient Kernel Discriminant
%   Analysis via Spectral Regression", Department of Computer Science
%   Technical Report No. 2888, University of Illinois at Urbana-Champaign
%   (UIUCDCS-R-2007-2888), August 2007.  
% 
%   6. Deng Cai, Xiaofei He, Jiawei Han, "Spectral Regression: A Unified
%   Subspace Learning Framework for Content-Based Image Retrieval", ACM
%   Multimedia 2007, Augsburg, Germany, Sep. 2007.
%
%   7. Deng Cai, Xiaofei He, Jiawei Han, "Spectral Regression for Efficient
%   Regularized Subspace Learning", IEEE International Conference on
%   Computer Vision (ICCV), Rio de Janeiro, Brazil, Oct. 2007. 
%
%   8. Deng Cai, Xiaofei He, Jiawei Han, "Spectral Regression: A Unified
%   Approach for Sparse Subspace Learning", Proc. 2007 Int. Conf. on Data
%   Mining (ICDM'07), Omaha, NE, Oct. 2007. 
%
%   9. Deng Cai, Xiaofei He, Jiawei Han, "Efficient Kernel Discriminant
%   Analysis via Spectral Regression", Proc. 2007 Int. Conf. on Data
%   Mining (ICDM'07), Omaha, NE, Oct. 2007. 
%
%   10. Deng Cai, Xiaofei He, Wei Vivian Zhang, Jiawei Han, "Regularized
%   Locality Preserving Indexing via Spectral Regression", Proc. 2007 ACM
%   Int. Conf. on Information and Knowledge Management (CIKM'07), Lisboa,
%   Portugal, Nov. 2007.
%
%
%   version 2.0 --Aug/2007 
%   version 1.0 --May/2006 
%
%   Written by Deng Cai (dengcai2 AT cs.uiuc.edu)
%


if ~exist('data','var')
    global data;
end

if ~exist('Responses','var')
    global Responses;
end



AppendConstant = 0;
if isfield(options,'AppendConstant') 
    AppendConstant = options.AppendConstant;
end

[nSmp,nFea] = size(data);
KernelWay = 0;
if nFea > nSmp
    KernelWay = 1;
end

nScale = min(nSmp,nFea);
if nScale < 3500
    if isfield(options,'LSQR') & options.LSQR
        options.LSQR = 0;
    end
end

nRepeat = 20;
if isfield(options,'nRepeat')
    nRepeat = options.nRepeat;
end

if ~isfield(options,'ReguAlpha')
    options.ReguAlpha = 0.1;
end


LassoCardi = 0;
if isfield(options,'LassoCardi') 
    LassoCardi = max(options.LassoCardi)+10;
end


tmp_T = cputime;

if isfield(options,'RemoveMean') & options.RemoveMean
    sampleMean = mean(data);
    data = (data - repmat(sampleMean,nSmp,1));
    AppendConstant = 0;
end

if AppendConstant
    data = [data ones(size(data,1),1)];
    [nSmp,nFea] = size(data);
end


if strcmpi(options.ReguType,'Lasso') 
    if nFea >= nSmp
        error('Please use RidgeLasso!');
    end
    options.ReguAlpha = 0;
    options.ReguType = 'RidgeLasso';
end


nVector = size(Responses,2);

itn = 0;
switch lower(options.ReguType)
    case {lower('Ridge')}
        if isfield(options,'LSQR') & options.LSQR
            [eigvector, istop, itn] = lsqr2(data, Responses, options.ReguAlpha, nRepeat);
        else
            if KernelWay
                ddata = data*data';
                ddata = full(ddata);
                ddata = single(ddata);

                for i=1:size(ddata,1)
                    ddata(i,i) = ddata(i,i) + options.ReguAlpha;
                end

                ddata = max(ddata,ddata');
                R = chol(ddata);
                eigvector = R\(R'\Responses);

                eigvector = double(eigvector);
                eigvector = data'*eigvector;
            else
                ddata = data'*data;
                ddata = full(ddata);
                ddata = single(ddata);

                for i=1:size(ddata,1)
                    ddata(i,i) = ddata(i,i) + options.ReguAlpha;
                end

                ddata = max(ddata,ddata');
                B = data'*Responses;

                R = chol(ddata);
                eigvector = R\(R'\B);
                eigvector = double(eigvector);
            end
        end
        eigvector = eigvector./repmat(max(1e-10,sum(eigvector.^2,1).^.5),size(eigvector,1),1);
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

        eigvector = {};
        if nFea < 5000
            global Gram;
            Gram = data'*data;
            Gram = max(Gram,Gram');
            
            if options.ReguAlpha > 0
                for i=1:size(Gram,1)
                    Gram(i,i) = Gram(i,i) + options.ReguAlpha;
                end
            end

            for i = ResponseRange
                disp(['Total projections: ',num2str(nVector),'; Projection now: ',num2str(i)]);
                eigvector_T = lars(data, Responses(:,i),'lasso', -LassoCardi,1,[],options.LassoCardi);
                if isfield(options,'tmpSave') & options.tmpSave
                    save(['eigvector_',num2str(i),'.mat'],'eigvector_T');
                else
                    eigvector{i} = eigvector_T;
                end
                clear eigvector_T;
            end
            clear global Gram;
        else
            if options.ReguAlpha > 0
                data = [data;sqrt(options.ReguAlpha)*speye(nFea)];
                Responses = [Responses;zeros(nFea,nVector)];
            end

            for i = ResponseRange
                disp(['Total projections: ',num2str(nVector),' Projection now: ',num2str(i)]);
                eigvector_T = lars(data, Responses(:,i),'lasso', -LassoCardi,0,[],options.LassoCardi);
                if isfield(options,'tmpSave') & options.tmpSave
                    save(['eigvector_',num2str(i),'.mat'],'eigvector_T');
                else
                    eigvector{i} = eigvector_T;
                end
                clear eigvector_T;
            end
            
            if options.ReguAlpha > 0
                data = data(1:nSmp,:);
                Responses = Responses(1:nSmp,:);
            end
            
        end
        
        if isfield(options, 'tmpExit') & options.tmpExit
            error('tmpExit');
        end
        
        if isfield(options,'tmpSave') & options.tmpSave
            for i = 1:nVector
                load(['eigvector_',num2str(i),'.mat'],'eigvector_T');
                eigvector{i} = eigvector_T;
            end
        end
    case {lower('Tensor')}
        ddata = data'*data;
        ddata = full(ddata);
        ddata = single(ddata);
        ddata = ddata + options.RegularizerOptions.ReguAlpha*options.RegularizerOptions.regularizerR;
        ddata = max(ddata,ddata');
        B = data'*Responses;

        R = chol(ddata);
        eigvector = R\(R'\B);
        eigvector = double(eigvector);

        eigvector = eigvector./repmat(max(1e-10,sum(eigvector.^2,1).^.5),size(eigvector,1),1);
    case {lower('Custom')}
        ddata = data'*data;
        ddata = full(ddata);
        ddata = single(ddata);
        ddata = ddata + options.RegularizerOptions.ReguAlpha*options.RegularizerOptions.regularizerR;
        ddata = max(ddata,ddata');
        B = data'*Responses;

        R = chol(ddata);
        eigvector = R\(R'\B);
        eigvector = double(eigvector);

        eigvector = eigvector./repmat(max(1e-10,sum(eigvector.^2,1).^.5),size(eigvector,1),1);
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
            tmpEigvec = eigvector_T(:,iMin(end))/norm(eigvector_T(:,iMin(end)));
            eigvector{cardidx} = [eigvector{cardidx} tmpEigvec];
        end
    end

    if nargout >= 4
        Y = cell(length(options.LassoCardi),1);
        for cardidx = 1:length(options.LassoCardi)
            Y{cardidx} = data*eigvector{cardidx};
        end
    end
else
    if nargout >= 4
        Y = data*eigvector;
    end
end


elapse = cputime - tmp_T;



