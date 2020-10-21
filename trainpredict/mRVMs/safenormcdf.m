function c = safenormcdf(x)
thresh=-10;
x(find(x<thresh))=thresh;
c=normcdf(x);