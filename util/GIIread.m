function g = GIIread(filename)

g = gifti(filename);
if ~isfield(g,'cdata')
    error('No cdata field found in GIFTI structure. Please make sure you read in a scalar map')
end