function x=choose_rnd(p)
% Chooses a random multinomial variable 
% with the probilities p
s=cumsum(p);
x=find(unifrnd(0,1)<s);
x=x(1);