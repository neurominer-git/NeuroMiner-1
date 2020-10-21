% wc_example
% 
% This script illustrates the usage of the Weibul-Cox functions using some
% synthetic data as an example.
%
% Copyright James Barrett 2014
%
% Version 1.1.4
% Date: 21 July 2014
% Contact: james.j.barrett@kcl.ac.uk

% Synthetic data generated elsewhere

data = [    1.8883    6.0000         0
    2.4348    6.0000         0
   -2.2381    0.8872         0
    2.4803    6.0000         0
    0.7942    1.4461         0
   -2.4148    1.1653    1.0000
   -1.3290    0.6275         0
    0.2813    2.9997    1.0000
    2.7450    6.0000         0
    2.7893    6.0000         0
   -2.0543    1.3937    1.0000
    2.8236    6.0000         0
    2.7430    6.0000         0
   -0.0877    3.2329    1.0000
    1.8017    4.4114    1.0000
   -2.1487    1.2023    1.0000
   -0.4694    1.5242         0
    2.4944    6.0000         0
    1.7532    6.0000         0
    2.7570    6.0000         0
    0.9344    4.1903    1.0000
   -2.7857    0.9338    1.0000
    2.0948    6.0000         0
    2.6040    6.0000         0
    1.0724    4.7664    1.0000];

X = data(:,1);
tau = data(:,2);
E = data(:,3);

%% Train model
model =  wc_train(X, tau, E)

%% Plot the inferred function
figure(1);clf;

% Use wc_predict to generate predictions with variance
x=-3:0.01:3;
p = zeros(length(x),1);     % predictions for each value of x
v = zeros(length(x),1);     % variance of these predictions
for i=1:length(x)
    [p(i), v(i)] = wc_predict(x(i),model);
end

% Plot these curves
L(1) = plot(x,p,'-k');
hold on
grid on
plot(x,p+sqrt(v),'--k')
plot(x,p-sqrt(v),'--k')

up = p+sqrt(v);    
Z=[x,fliplr(x)];          
W=[p',fliplr(up')];              
fill(Z,W,[0.8 0.8 0.8]);                  

lp = p-sqrt(v);    
Z=[x,fliplr(x)];          
W=[p',fliplr(lp')];              
fill(Z,W,[0.8 0.8 0.8]); 

% Plot the event times
data = [X,tau,E];
data = sortrows(data,3);
X0 = data(1:length(find(model.E==0)),1);
tau0 = data(1:length(find(model.E==0)),2);
L(2) = plot(model.X,model.t,'ko','MarkerFaceColor','k');
L(3) = plot(X0,tau0,'ko','MarkerFaceColor','w');
xlabel('Covariate x')
ylabel('Time to event t')
title('Example of Cox PH model with Weibull base hazard')
legend(L,'Inferred function','Non-censored events','Censored events','location','northwest')
set(gca(),'Xlim',[-3,3])

%% Generate hazard rate and survival funciton for xstar=1
t = 0:0.1:10;
xstar = 1;
h = wc_hazard(t,xstar,model);
S = wc_survival(t,xstar,model);

figure(2);clf;
plot(t,S,'-k')
grid on
xlabel('Time t')
ylabel('Survival probability')
title('Survival function for xstar = 1')

figure(3);clf;
plot(t,h,'-k')
set(gca(),'Xlim',[0,5])
grid on
xlabel('Time t')
ylabel('Hazard rate')
title('Hazard rate for xstar = 1')