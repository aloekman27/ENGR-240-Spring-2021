% This script analyses the use of 3 ODE solver which are 4th Order Runge-Kutta,
% ode45 and ode23s with Tolerances at 1e-8.

%define the constants
y0 = [1.0 2.2 0]; %Initial conditions A0,B0,C0
tspan = [0 6]; %time span
k1 = 0.25; k2 = 0.10; k3 = 125; %k values 
h = [0.015 0.01 0.005]; %time step values, h


%define the anonymous function
dydt = @(t,y,k1,k2,k3) [(-k1.*y(1).*y(2))+(k2.*y(3));(-k1.*y(1).*y(2))+(k2.*y(3))-(k3.*y(2).*y(3));(k1.*y(1).*y(2))-(k2.*y(3))-(k3.*y(2).*y(3))];

%1st output uses 4th Order Runge-Kutta with h = 0.015
[t1,ABC1] = odeRK4sys(dydt,tspan,y0,h(1),k1,k2,k3);

%2nd output uses 4th Order Runge-Kutta with h = 0.01
[t2,ABC2] = odeRK4sys(dydt,tspan,y0,h(2),k1,k2,k3);

%3rd output uses 4th Order Runge-Kutta with h = 0.005
[t3,ABC3] = odeRK4sys(dydt,tspan,y0,h(3),k1,k2,k3);

%set options to be used in ode45 and ode23s
opts = odeset('RelTol',1e-8,'AbsTol',1e-8);

%4th output uses MATLAB's ode45 function 
[t4,ABC4] = ode45(dydt,tspan,y0,opts,k1,k2,k3);

%5th output uses MATLAB's ode23s function 
[t5,ABC5] = ode23s(dydt,tspan,y0,opts,k1,k2,k3);

%find each time step of every tvalues from the evaluation of the function
T45 = diff(t4);
T23s = diff(t5);

%find the mean t values which is the average time step used in the function evaluation
Tmean45 = mean(T45);
Tmean23s = mean(T23s);