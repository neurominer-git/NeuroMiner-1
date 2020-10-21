function [D, mD] = nk_ComputeMnMdPairWiseDiff(G, func, act)

switch func
    case 'mn'
        func = @nm_nanmean;
    case 'md'
        func = @nm_nanmedian;
end
[m,n] = size(G);

switch act
    case 'meandiff'
        D = zeros(m,n);
        for i=1:n
            vec = true(1,n); vec(i)=false;
            D(:,i) = func(G(:,i)-G(:,vec),2);
        end
        mD = func(D);
    case 'alldiff'
        D = zeros(m,n-1,n);
        mD = zeros(1,n-1,n);
        for i=1:n
            vec = true(1,n); vec(i)=false;
            D(:,:,i) = G(:,i)-G(:,vec);
            mD(:,:,i) = func(D(:,:,i));
        end
        
end