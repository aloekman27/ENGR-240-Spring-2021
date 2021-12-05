%% Seasonal Spread of Infectious Disease ODE and Parameter Analysis
%	This script m-file solves a system of differential equations that
%	is proposed as a model of the seasonal spread of infectious disease in
%	a population. A function B(t) is also used in the differential equation
%	to describe the seasonally dependent infection rate modeled as the
%	periodic function.
%
% Time span is 48 months
% Initial Conditions:
%       S(0) = 1    S = fractions of susceptible individuals in population
%       I(0) = 1e-6 I = fractions of infected individuals in population
%       R(0) = 0    R = fractions of recovered individuals in population
% Season length, T = 12
% Initial infection rate b0 = 3.0
% Birth death rate mu = 0.1
% Recovery rate gamma = 1.5
% 
% Outputs are in the form of t and y where the 1st column of y is S, 2nd
% column of y is I, and 3rd column of y is R


%Defining all the parameter and initial conditions for the system of ODE
y0 = [1 1e-6 0]; %initial conditions S(0), I(0), and R(0)
T = 12; %season length
b_0 = 3.0; %initial infection rate
mu = [0.1 0.14 0.25 0.37 0.50]; %birth date rate
gamma = [1.5 0.001 0.8 1.8 3]; %recovery rate
tspan = [0 48]; %time span in months
t0 = tspan(1); tf = tspan(end); %The start and final time interval
opts = odeset('RelTol',1e-8,'AbsTol',1e-8); %ODE solver options

%Anonymous function for the seasonally dependent infection rate
b_t = @(t) b_0*(1+sin((2*pi*t)/T));

%Anonymous functions of the seasonal spread of infectious disease in a
%population where y(1) is S (susceptible), y(2) is I (infected), and y(3)
%is R (recovered) individuals
dydt = @(t,y,b_t,mu,gamma) [mu.*(1-y(1))-(b_t(t).*y(2).*y(1));b_t(t).*y(2).*y(1)-(gamma+mu).*y(2);gamma.*y(2)-mu.*y(3)];

%% SOLVING THE SYSTEM OF ODES FOR t = 0 to 48 MONTHS
% This section solves the system of ODEs using different methods which are
% ode45, ode23s, and odeRK4sys

%Solving the system of ODE using a built-in function ode45
tic
[t1,y1] = ode45(@(t,y)dydt(t,y,b_t,mu(1),gamma(1)),tspan,y0,opts);
toc
%Finding the average step size of the time interval taken by the built-in
%function ode45
Tmean45 = mean(diff(t1));

%Plotting the values of S, I, and R over the time interval that are solved
%using a built-in function ode45 on a single graph
figure(1)
plot(t1,y1)
title('The graph of S,I,R using ode45')
xlabel('time (months)')
ylabel('fraction of individuals in population')
legend('S','I','R')
grid

%Solving the system of ODE using a built-in function ode23s
tic
[t2,y2] = ode23s(@(t,y)dydt(t,y,b_t,mu(1),gamma(1)),tspan,y0,opts);
toc
%Finding the average step size of the time interval taken by the built-in
%function ode23s
Tmean23s = mean(diff(t2));

%Plotting the values of S, I, and R over the time interval that are solved
%using a built-in function ode23s on a single graph
figure(2)
plot(t2,y2)
title('The graph of S,I,R using ode23s')
xlabel('time (months)')
ylabel('fraction of individuals in population')
legend('S','I','R')
grid

%Proof that the susceptible (S), infected (I), and recovered (R)
%individuals are fractions to a population resulting in a total value of 1
S = y1(:,1);
I = y1(:,2);
R = y1(:,3);
test1 = S + I + R;

%Proof that the total of the derivative of susceptible (S), infected (I),
%and recovered (R) individuals is 0
dSdt = gradient(S,Tmean45);
dIdt = gradient(I,Tmean45);
dRdt = gradient(R,Tmean45);
test2 = dSdt + dIdt + dRdt;

%Plotting the graph of both tests to proof that both values remain constant
%over the time interval
figure(8)
plot(t1,test1,t1,test2)
title('showing that S + I + R = 1 and dSdt + dIdt + dRdt = 1')
legend('S + I + R','dSdt + dIdt + dRdt')

%Solving the system of ODE using a 4th order Runge-Kutta method
tic
[t3,y3] = odeRK4sys(@(t,y)dydt(t,y,b_t,mu(1),gamma(1)),tspan,y0,0.01);
toc

%Plotting the values of S, I, and R over the time interval that are solved
%using 4th order Runge-Kutta method on a single graph
figure(3)
plot(t3,y3)
title('The graph of S,I,R using 4th order Runge-Kutta method')
xlabel('time (months)')
ylabel('fraction of individuals in population')
legend('S','I','R')
grid

%% Paraneter Sensitivity Analysis 1 (Varying Gamma)

%Solving the ODE system for different recovery rate (gamma) values using a 
%built-in function ode45
%Gamma = 0.001
[t12,y12] = ode45(@(t,y)dydt(t,y,b_t,mu(1),gamma(2)),tspan,y0,opts);
%Gamma = 0.8
[t13,y13] = ode45(@(t,y)dydt(t,y,b_t,mu(1),gamma(3)),tspan,y0,opts);
%Gamma = 1.8
[t14,y14] = ode45(@(t,y)dydt(t,y,b_t,mu(1),gamma(4)),tspan,y0,opts);
%Gamma = 3
[t15,y15] = ode45(@(t,y)dydt(t,y,b_t,mu(1),gamma(5)),tspan,y0,opts);

%Finding the I max for each of the gamma value
%I max value and the corresponding time for gamma = 1.5
[maxvaly1,idx1] = max(y1(:,2));
maxvalx1 = t1(idx1);
%I max value and the corresponding time for gamma = 0.001
[maxvaly2,idx2] = max(y12(:,2));
maxvalx2 = t12(idx2);
%I max value and the corresponding time for gamma = 0.8
[maxvaly3,idx3] = max(y13(:,2));
maxvalx3 = t13(idx3);
%I max value and the corresponding time for gamma = 1.8
[maxvaly4,idx4] = max(y14(:,2));
maxvalx4 = t14(idx4);
%I max value and the corresponding time for gamma = 3
[maxvaly5,idx5] = max(y15(:,2));
maxvalx5 = t15(idx5);

%Plotting all the values of I and t with different gamma values that are
%obtained from ode45 function and indicating the I max value in each graph
figure(4)
hold on
plot(t1,y1(:,2),'LineWidth',1)
plot(t12,y12(:,2),'--')
plot(t13,y13(:,2),'--')
plot(t14,y14(:,2),'--')
plot(t15,y15(:,2),'--')
plot(maxvalx1,maxvaly1,'rd','LineWidth',1)
plot(maxvalx2,maxvaly2,'rd','LineWidth',1)
plot(maxvalx3,maxvaly3,'rd','LineWidth',1)
plot(maxvalx4,maxvaly4,'rd','LineWidth',1)
plot(maxvalx5,maxvaly5,'rd','LineWidth',1)
title('Gamma parameter sensitivity analysis using ode45')
xlabel('Time (Months)')
ylabel('Infected Population (I)')
legend ('gamma = 1.5','gamma = 0.001','gamma = 0.8','gamma = 1.8','gamma = 3')
hold off

%Solving the ODE system for different gamma values using a built-in 
%function ode23s
%Gamma = 0.001
[t112,y112] = ode23s(@(t,y)dydt(t,y,b_t,mu(1),gamma(2)),tspan,y0,opts);
%Gamma = 0.8
[t113,y113] = ode23s(@(t,y)dydt(t,y,b_t,mu(1),gamma(3)),tspan,y0,opts);
%Gamma = 1.8
[t114,y114] = ode23s(@(t,y)dydt(t,y,b_t,mu(1),gamma(4)),tspan,y0,opts);
%Gamma = 3
[t115,y115] = ode23s(@(t,y)dydt(t,y,b_t,mu(1),gamma(5)),tspan,y0,opts);

%Plotting all the values of I and t with different gamma values that are
%obtained from ode23s function
figure(5)
hold on
plot(t2,y2(:,2),'LineWidth',1)
plot(t112,y112(:,2),'--')
plot(t113,y113(:,2),'--')
plot(t114,y114(:,2),'--')
plot(t115,y115(:,2),'--')
title('Gamma parameter sensitivity analysis using ode23s')
xlabel('Time (Months)')
ylabel('Infected Population (I)')
legend ('gamma = 1.5','gamma = 0.001','gamma = 0.8','gamma = 1.8','gamma = 3')
hold off

%% Paraneter Sensitivity Analysis 1 (Varying mu)

%Solving the ODE system for different birth date rate (mu) values using a 
%built-in function ode45
%mu = 0.14
[t22,y22] = ode45(@(t,y)dydt(t,y,b_t,mu(2),gamma(1)),tspan,y0,opts);
%mu = 0.25
[t23,y23] = ode45(@(t,y)dydt(t,y,b_t,mu(3),gamma(1)),tspan,y0,opts);
%mu = 0.46
[t24,y24] = ode45(@(t,y)dydt(t,y,b_t,mu(4),gamma(1)),tspan,y0,opts);
%mu = 0.50
[t25,y25] = ode45(@(t,y)dydt(t,y,b_t,mu(5),gamma(1)),tspan,y0,opts);

%Finding the I max for each of the mu value
%I max value and the corresponding time for mu = 0.14
[maxvaly6,idx6] = max(y22(:,2));
maxvalx6 = t22(idx6);
%I max value and the corresponding time for mu = 0.25
[maxvaly7,idx7] = max(y23(:,2));
maxvalx7 = t23(idx7);
%I max value and the corresponding time for mu = 0.46
[maxvaly8,idx8] = max(y24(:,2));
maxvalx8 = t24(idx8);
%I max value and the corresponding time for mu = 0.50
[maxvaly9,idx9] = max(y25(:,2));
maxvalx9 = t25(idx9);

%Plotting all the values of I and t with different mu values that are
%obtained from ode45 function and indicating the I max value in each graph
figure(6)
hold on
plot(t1,y1(:,2),'LineWidth',1)
plot(t22,y22(:,2),'--')
plot(t23,y23(:,2),'--')
plot(t24,y24(:,2),'--')
plot(t25,y25(:,2),'--')
plot(maxvalx1,maxvaly1,'rd','LineWidth',1)
plot(maxvalx6,maxvaly6,'rd','LineWidth',1)
plot(maxvalx7,maxvaly7,'rd','LineWidth',1)
plot(maxvalx8,maxvaly8,'rd','LineWidth',1)
plot(maxvalx9,maxvaly9,'rd','LineWidth',1)
title('Mu parameter sensitivity analysis using ode45')
xlabel('Time (Months)')
ylabel('Infected Population (I)')
legend('mu = 0.1','mu = 0.14','mu = 0.25','mu = 0.46','mu = 0.50')
hold off

%Solving the ODE system for different mu values using a built-in 
%function ode23s
%mu = 0.14
[t212,y212] = ode23s(@(t,y)dydt(t,y,b_t,mu(2),gamma(1)),tspan,y0,opts);
%mu = 0.25
[t213,y213] = ode23s(@(t,y)dydt(t,y,b_t,mu(3),gamma(1)),tspan,y0,opts);
%mu = 0.46
[t214,y214] = ode23s(@(t,y)dydt(t,y,b_t,mu(4),gamma(1)),tspan,y0,opts);
%mu = 0.50
[t215,y215] = ode23s(@(t,y)dydt(t,y,b_t,mu(5),gamma(1)),tspan,y0,opts);

%Plotting all the values of I and t with different mu values that are
%obtained from ode23s function
figure(7)
hold on
plot(t2,y2(:,2),'LineWidth',1)
plot(t212,y212(:,2),'--')
plot(t213,y213(:,2),'--')
plot(t214,y214(:,2),'--')
plot(t215,y215(:,2),'--')
title('Mu parameter sensitivity analysis using ode23s')
xlabel('Time (Months)')
ylabel('Infected Population (I)')
legend('mu = 0.1','mu = 0.14','mu = 0.25','mu = 0.46','mu = 0.50')
hold off
