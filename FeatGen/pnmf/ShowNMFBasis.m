function ShowNMFBasis(W, w, h, scale, dotranspose)
[m, r] = size(W);
rsq = sqrt(r);

if ~exist('scale', 'var') || isempty(scale)
    scale = 2;
end

if ~exist('dotranspose', 'var') || isempty(dotranspose)
    dotranspose = false;
end

w1 = w * scale;
h1 = h * scale;
padding = 3;

titleheight = 0;
allwidth = w1 * rsq + padding * (rsq - 1);
allheight = h1 * rsq + padding * (rsq - 1) + titleheight;

figure('Position', [100,100,allwidth,allheight]);
for i=1:rsq
    for j=1:rsq
        t = (i-1) * rsq + j;
        subplot('Position', [(w1+padding)*(j-1)/allwidth,(h1+padding)*(rsq-i)/allheight,w1/allwidth,h1/allheight]);
        colormap gray;
        if dotranspose
            imagesc(reshape(W(:,t), [h,w]));
        else
            imagesc(reshape(W(:,t), [w,h])');
        end

        axis off;
    end
end
