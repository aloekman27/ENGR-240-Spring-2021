%Arrhenius Equation 2

% Use Non-linear Regression Method (Fminsearch) and General Linear Least
% Squares to determine values of A and E that achieves the least squares
% fit for the modified Arrhenius model to the data
%                    k = A.*(T_a.^b).*exp(-E./(R.*T_a))

load O_H2_rxn.dat %load the data from .txt file 
tic
%Set up known values used in the equation
R = 8.314; %J/(mol*K)
T_a = O_H2_rxn(:,1)'; %Kelvin, K
k = O_H2_rxn(:,2)'; %s^-1

%% Nonlinear Regression (fminsearch)
% % %  This part of the script m-file uses the Nonlinear Regression
% (fminsearch) method to curve fit the Arrhenius Equation
%                           k = A.*exp(-E./(R.*T_a))
%

%Define the equation for the model
k_model = @(T_a,A,b,E) A.*(T_a.^b).*exp(-E./(R.*T_a));

%Find the squared to norm of the residual to use fminsearch
k_res = @(a,T_ai,ki) sum((ki-k_model(T_ai,a(1),a(2),a(3))).^2);
options = optimset('TolX',1e-6,'TolFun',1e-6);
a = fminsearch(k_res,[40 35 45],options,T_a,k);
A_nonl = a(1)
b_nonl = a(2)
E_nonl = a(3)

%Calculate the quality of the fit
St = sum((k-mean(k)).^2);
Sr_nonl = sum((k-k_model(T_a,A_nonl,b_nonl,E_nonl)).^2);
r_squared_nonl = 1-Sr_nonl/St
StandardError_nonl = sqrt(Sr_nonl/(length(k)-length(a)))
toc

%% General Linear Least Squares
%  This part of the script m-file uses the General Linear Least Squares
%  method to curve fit the Arrhenius Equation
%                           k = A.*exp(-E./(R.*T_a))
%
tic
%define linearized y
y_GLLS = log(k);

%basis functions
Z0 = @(x) x.^0;
Z1 = @(x) log(x);
Z2 = @(x) (1./(x.*R));

%put basis functions with input to a matrix
Z = [Z0(T_a)',Z1(T_a)',Z2(T_a)'];

%calculate coefficients
abc = Z\y_GLLS';

%reverse linearization to calculate constants
A_GLLS = exp(abc(1))
b_GLLS = abc(2)
E_GLLS = -abc(3)

%preallocate values to calculate standard error
d = length(abc);
n = length(k);

%calculate the quality of the fit
Sr_GLLS = sum((k-k_model(T_a,A_GLLS,b_GLLS,E_GLLS)).^2);%sum of the squared residuals
rSquared = 1-Sr_GLLS/St
StandardError_GLLS = sqrt(Sr_GLLS/(n-d))

%Preallocate variables that will help to plot the graph
T_a_plot = linspace(min(T_a),max(T_a),100);
yplot = log(A_GLLS)*Z0(T_a_plot)+b_GLLS*Z1(T_a_plot)-E_GLLS*Z2(T_a_plot);
kplot = exp((yplot));
toc

%% Plot the Graph
% This is the figure for the curve fitting method we have done above. 

%Plotting the graph
T_am = min(T_a):0.1:max(T_a);
figure(1)
plot(T_a,k,'ko',T_am,k_model(T_am,A_nonl,b_nonl,E_nonl),'m--',T_a_plot,kplot,'b--')
title('Temperature Dependence of the Reaction Rate')
xlabel('Temperature (K)')
ylabel('Reaction Rate (s^{-1})')
legend('data','fminsearch','General Linear Least Squares')
