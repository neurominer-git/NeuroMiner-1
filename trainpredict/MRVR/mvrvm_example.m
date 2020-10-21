%Demonstrates the mvrvm algorithm on a toy example. Part of the code was
%derived from Tipping's matlab code for top down rvm for sing variate
%output.


%  @file mvrvm_example.m
%  				Author Arasanathan Thayananthan ( at315@cam.ac.uk)
%               (c) Copyright University of Cambridge
%  
%     This library is free software; you can redistribute it and/or
%     modify it under the terms of the GNU Lesser General Public
%     License as published by the Free Software Foundation; either
%     version 2 of the License, or (at your option) any later version.
%  
%     This library is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%     Lesser General Public License for more details.
%  
%     You should have received a copy of the GNU Lesser General Public
%     License along with this library; if not, write to the Free Software
%     Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 
 



function MVRVM_example



close all;
 randn('state',5)

 N		= 101;
 noise		= [0.1 0.1];
 kernel_	= '+Gauss';
 width		= 3.0;
 maxIts	= 10;

% Generate the sinc data
 
X	= 10*[-1:2/(N-1):1]';

y1	= sin(abs(X))./abs(X);
y1(find(X==0))=1;
% Generate linear data
y2=0.1*X+1.5;


%Multivariate output
y=[y1 y2]

%Introduce some noise
t	= y + repmat(noise,N,1).*randn(N,2);


% Plot the data and function
% 
figure(1)
hold on;
h_y=plot(X,y(:,1),'r--','LineWidth',2);
plot(X,y(:,2),'r--','LineWidth',2);
plot(X,t(:,1),'black.','MarkerSize',16)
plot(X,t(:,2),'black.','MarkerSize',16)

box = [-10.1 10.1 1.1*[min(min(t)) max(max(t))]];
axis(box)
drawnow


PHI	= sbl_kernelFunction(X,X,kernel_,width);


[weights, used, alpha, beta] = mvrvm(PHI,t,maxIts);

%need to check if bias is used or not
 
if kernel_(1)=='+'	
  used	= used - 1;
  if used(1)~=0
    kernel_(1)	= [];
  else
    used(1)	= [];
  end
end

 

%
% Plot the results
PHI	= sbl_kernelFunction(X,X(used,:),kernel_,width);
y_rvm	= PHI*weights; 

h_yrvm = plot(X,y_rvm(:,1),'b-','LineWidth',2);
h_rv	= plot(X(used),t(used,1),'go','LineWidth',2,'MarkerSize',10);
h_yrvm = plot(X,y_rvm(:,2),'b-','LineWidth',2);
h_rv	= plot(X(used),t(used,2),'go','LineWidth',2,'MarkerSize',10);
legend([h_y h_yrvm h_rv],'target functions','RVM predictor','RVs')
hold off
%title('RVM for Multivariate Traget Functions','FontSize',14)
%
% Output some info
% 
fprintf('\nRVM regression test error (RMS): %g\n', ...
	sqrt(mean((y-y_rvm).^2)))
fprintf('\testimated noise level: %.4f (true: %.4f)\n', ...
	sqrt(1./beta), noise)