function IO = ReadSurf(IO)

cnt = 1; IO.Y=[];
for j=1:numel(IO.V)
    for i=1:IO.n_subjects(j)
        filename = deblank(IO.PP(cnt,:));
        [~,y] = SurfaceReader(filename);
        fprintf('\nsample %g, file: %s',j, filename)
        if size(y,1)>1, y=y';end
        IO.Y = [IO.Y; y];
        cnt=cnt+1;
    end
end