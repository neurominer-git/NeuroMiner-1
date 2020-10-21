function y=fcneval(fcn,x)
% Evaluates a function, while dealing with possible complex strings 
% like sqrt(x.^2)
% x is the placeholder for the data 
if (any(fcn=='('))
    y=eval(fcn);        
else 
    y=feval(fcn,x);
end;
