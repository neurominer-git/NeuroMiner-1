function G = generateDAG(T)

[mT,nT]=size(T); 
G = zeros(mT,mT-1);
G(1,1)=1;
i=2; 
oT=T; 
while i<mT
    cnt = 1;
    d = G(:,i-1);
    dvec = find(d);
    for j=1:numel(dvec)
        jT = oT(:,d(dvec(j)));
        tT = oT(:,d(dvec(j))+1:nT);
        Ip = find( tT( jT ==  1, :) == 1 );
        In = find( tT (jT == -1, :) == 1 );
        [~,Fp] = ismember( tT(:,Ip(1))', oT', 'rows' );
        [~,Fn] = ismember( tT(:,In(1))', oT', 'rows' );
        G(cnt,i) = Fp; cnt = cnt+1;
        G(cnt,i) = Fn; cnt = cnt+1;
    end
    i=i+1;
end 