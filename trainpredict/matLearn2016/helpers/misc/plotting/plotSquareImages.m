function [] = plotSquareImages(images)
% Uses method subtightplot from Matlab community, please see licence file 
% for usage rights

[M, L]=size(images);
squareSide = sqrt(L);

if squareSide^2 ~= L
    error('Images must be square')
end

rowLen = ceil(sqrt(M));
colLen = ceil(sqrt(M));
figure;
thinMargin = [0.005 0.005];
subplot = @(m,n,p) subtightplot (m, n, p, thinMargin, thinMargin, ...
                                 thinMargin);
for i = 1:M
            im = reshape(images(i,:),[squareSide, squareSide]);   
            h = subplot(rowLen,colLen,i);
            imagesc(im);
            set(h, 'Visible', 'off');
end

end