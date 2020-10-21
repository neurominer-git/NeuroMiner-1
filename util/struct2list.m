function L=struct2list(M)
% function L=struct2list(M)
% Takes the field of a strcuture and converts them into a argument list    
S=fieldnames(M);
    for i=1:length(S)
        L{i*2-1}=S{i};
        L{i*2}=getfield(M,S{i});
    end;
    