function [d, v] = SurfaceReader(filename)

[~,~,e] = fileparts(deblank(filename));
switch e
    case {'.mgh','.mgz'}
        d = MRIread(filename); v = d.vol;
    case '.gii'
        d = GIIread(filename); v = d.cdata'; 
    otherwise
        error('Do not know how to process surface files with %s extension',e)
end