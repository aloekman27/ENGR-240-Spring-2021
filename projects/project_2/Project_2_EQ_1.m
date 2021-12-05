%Arrhenius Equation 1

% Use curve fitting methods, Linear Regression and General Linear Least
% Squares, to determine values of A and E that achieves the least squares
% fit for the unmodified Arrhenius model to the data

clear all 
clc
load O_H2_rxn.dat %load the data from .txt file 

%% Linear Regression
% This part of the script m-file uses linear regression method to curve fit
% the Arrhenius Equation 
%                           k = A.*exp(-E./(R.*T_a))
%
tic
%Set up known values used in the equation
R = 8.314; %J/(mol*K)
T_a = O_H2_rxn(:,1)'; %Kelvin, K
k = O_H2_rxn(:,2)'; %s^-1

%Linear Regression of transformed model
p = polyfit(1./(R.*T_a),log(k),1);
A1 = exp(p(2))
E1 = -p(1)
kmod_lin = @(T_a) A1.*exp(-E1./(R.*T_a)); 

%calculate the quality of the fit
Sr_lin = sum((k-kmod_lin(T_a)).^2);
St_lin = sum((k-mean(T_a)).^2);
rSquared_lin = 1-Sr_lin/St_lin
StandardError_lin = sqrt(Sr_lin/(length(k)-2))


%preallocate values for the graph later
T_a_m = min(T_a):0.1:max(T_a);
toc
%% General Linear Least Square
%  This part of the script m-file uses the General Linear Least Squares
%  method to curve fit the Arrhenius Equation
%                           k = A.*exp(-E./(R.*T_a))
%
tic
%define the linearized equation of y
y_GLLS = log(k);

%basis functions
Z0 = @(x) x.^0;
Z1 = @(x) -(1./(R.*x));
Z = [Z0(T_a)',Z1(T_a)'];

%do left division to find the value of A and E, for GLLD
A = Z\y_GLLS';
A2 = exp(A(1))
E2 = A(2)


d = length(A); %degree of freedom
n = length(k); % number of data

%set up the k model for the equation
k_model = @(T_a) A2.*exp(-E2./(R.*T_a));

%Quality of the fit for General Linear Least squares approach
St_GLLS = sum((k-mean(k)).^2); % the measure of St in GLLS
Sr_GLLS = sum((k-k_model(T_a)).^2);%sum of the squared residuals, GLLS
rSquared_GLLS = 1-Sr_GLLS/St_GLLS

%Standard error of the fit
standardError_GLLS = sqrt(Sr_GLLS/n-d)

%unreverse the lienarization and plot the data
T_a_plot = linspace(min(T_a),max(T_a),100);
yplot = log(A2)*Z0(T_a_plot)+E2*Z1(T_a_plot);
kplot = exp(yplot);
toc

%% NonLinear Regression (fminsearch)
% %  This part of the script m-file uses the Nonlinear Regression
% (fminsearch) method to curve fit the Arrhenius Equation
%                           k = A.*exp(-E./(R.*T_a))
%
tic
%define the equation for the model
kmodel = @(T_a,A,E) A.*exp(-E./(R.*T_a));

%find the squared to norm of the residual to use fminsearch
k_res = @(a,T_ai,ki) sum((ki-kmodel(T_ai,a(1),a(2))).^2);
options = optimset('TolX',1e-2,'TolFun',1e-2);
a = fminsearch(k_res,[40,40],options,T_a,k);
A_nonl = a(1)
E_nonl = a(2)

%calculate the quality of the fit (fminsearch)
St_nonl = sum((k-mean(k)).^2);
Sr_nonl = sum((k-kmodel(T_a,A_nonl,E_nonl)).^2);
r_squared_nonl = 1-Sr_nonl/St_nonl
StandardError_nonl = sqrt(Sr_nonl/(length(k)-2))


T_am = min(T_a):0.1:max(T_a);
toc


%% Plot the Graph
% This is the figure for the curve fitting method we have done above. 

figure(1)
plot(T_a,k,'kd',T_a_m,kmod_lin(T_a_m),'b--',T_a_plot,kplot,'g--',T_am,kmodel(T_am,A_nonl,E_nonl),'r--')
title('Temperature Dependence of the Reaction Rate')
xlabel('Temperature (K)')
ylabel('Reaction Rate (s^{-1})')
legend('data','linear regression','General Linear Least Square','fminsearch')


