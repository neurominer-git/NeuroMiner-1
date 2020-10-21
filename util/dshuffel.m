function dshuffel(infile,outfile) 
D=dload(infile); 
F=fieldnames(D);
N=size(D.(F{1}),1);
D=getrow(D,sample_wor([1:N],N));
dsave(outfile,D);
