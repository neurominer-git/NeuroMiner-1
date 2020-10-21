function [meanInner,medianInner,ARI,sortedBasisNum] = evaluateReproducibility(pathDir1,pathDir2)

% function that evaluates the reproducibility between the results obtained
% for two different splits
% The scripts assumes that experiments using the same range of components
% have been performed for both splits

% outputs
% meanInner : mean value of the inner product between matched components
% medianInner : median value of the inner product between matched components
% ARI : adjusted Rand Index evaluated by deriving hard clusters from the
% estimated components
% sortedBasisNum : range of values for which components were estimated

listing = dir(pathDir1);
listing=listing(3:end) ;
hh =cellfun(@(x) ( (strfind(x,'NumBases')==1)  ),{listing(:).name},'UniformOutput',false) ;
listing=listing(cellfun(@(x) (~isempty(x)&&(x==1)),hh));
numDifBases=numel(listing);

% sort them in ascending order
basisNum = zeros(1,numDifBases) ;
for i=1:numDifBases
    basisNum(i) = str2double(listing(i).name(9:end));
end
[~,idx]=sort(basisNum) ;
sortedBasisNum=basisNum(idx) ;

ARI=zeros(numDifBases,1);

for exp=1:numDifBases
    disp([ num2str(exp) '/' num2str(numDifBases)])
   
    resSplit1 = load([pathDir1 '/NumBases' num2str(sortedBasisNum(exp)) '/OPNMF/ResultsExtractBases.mat']) ;
    resSplit2 = load([pathDir2 '/NumBases' num2str(sortedBasisNum(exp)) '/OPNMF/ResultsExtractBases.mat']) ;
    
    % normalize to unit norm
    blen1 = sqrt(sum((resSplit1.B).^2)) ;
    blen2 = sqrt(sum((resSplit2.B).^2)) ;    
    
    if any(blen1==0)
        blen1(blen1==0) = 1;
    end
    B1 = bsxfun(@times,resSplit1.B,1./blen1) ;
   
    if any(blen2==0)
        blen2(blen2==0) = 1;
    end
    B2 = bsxfun(@times,resSplit2.B,1./blen2) ;
    

    % calculate inner products
    inner_product = B1'*B2 ;
    
    % take a distance
    dist = 2*(1 - inner_product) ;
    
    % find correspondences
    [Matching,~] = Hungarian(dist);
    [~,idx_hug1]=max(Matching,[],2);
    
    % overlap - hungarian
    overlap{exp} = zeros(length(blen1),1) ;
    for b=1:length(blen1)
        overlap{exp}(b) = inner_product(b,idx_hug1(b));
    end
    
    % overlap with best
    overlap_best{exp} = max(inner_product,[],2) ;
    
    % also evaluate overlap based on adjusted Rand Index    
    rowLen1 = sum(B1,2) ;
    rowLen2 = sum(B2,2) ;
    
    if any(rowLen1==0)
        rowLen1(rowLen1==0) = 1 ;
    end
    if any(rowLen2==0)
        rowLen2(rowLen2==0) = 1 ;
    end
    BB1 = bsxfun(@times,(B1'),1./(rowLen1')); BB1=BB1';
    BB2 = bsxfun(@times,(B2'),1./(rowLen2')); BB2=BB2';
    
    [~,clustering1] = max(BB1,[],2);
    [~,clustering2] = max(BB2,[],2);
    ARI(exp) = clustering_adjustedRand_fast(clustering1,clustering2);
    
end

meanInner=cellfun(@(x) mean(x),overlap,'UniformOutput', false);
medianInner=cellfun(@(x) median(x),overlap,'UniformOutput', false);

%stdInner=cellfun(@(x) std(x),overlap,'UniformOutput', false);
figure;plot(sortedBasisNum,cell2mat(meanInner),'b','LineWidth',2)
xlabel('Number of components','fontsize',12)
ylabel({'Split-sample reproducibility';'(mean inner product)'},'fontsize',12)
set(gca,'fontsize',12)


figure;plot(sortedBasisNum,cell2mat(medianInner),'b','LineWidth',2)
xlabel('Number of components','fontsize',12)
ylabel({'Split-sample reproducibility';'(median inner product)'},'fontsize',12)
set(gca,'fontsize',12)

figure;plot(sortedBasisNum,ARI,'b','LineWidth',2)
xlabel('Number of components','fontsize',12)
ylabel({'Split-sample reproducibility' ;'(Adjusted Rand Index)'},'fontsize',12)
set(gca,'fontsize',12)