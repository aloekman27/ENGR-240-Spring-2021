%example function call for first test
T = 30200; %N
C_D = 0.088; %N*s^2/m^2
BR = 30; %kg/s
tspan = [0 20];
[time, altitude, velocity] = PA8_rocket(T,C_D,BR,tspan,0.5);

function [time, altitude, velocity] = PA8_rocket(thrust, drag_coeff, fuel_burn_rate, tspan, h)
% [time, altitude, velocity] = PA8_rocket(thrust, drag_coeff, fuel_burn_rate, tspan, h)
%
% Use midpoint method to solve a system of ODEs and model the kinematics of a rocket launch
%
% Input:
%   thrust         : A scalar input for the thrust (Newtons)
%   drag_coeff     : A scalar input for the drag coefficient (kg/m)
%   fuel_burn_rate : A scalar input for the fuel burn rate (kg/s)
%   tspan          : two-elemet vector for the initial and final time for the ODE
%   h              : a scalar input for the step size 
%
% Output:
%   time           : a column vecotr of time values generated from the tspan & h inputs 
%   altitude       : a column vector of altitude corresponding to the time values 
%   velocity       : a column vector of altitude corresponding to the time values

g = 9.81; %acceleration due to gravity constant (m/s^2)

%define inputs to a better variable
T = thrust; %Newtons
cd = drag_coeff; %kg/m
fbr = fuel_burn_rate;%kg/s

% define an anonymous function for the ODE which is to be solved 
dydt=@(t,y)[y(2);(T - (1200-fbr*t)*g-cd*(y(2)).^2)/(1200-fbr*t)];

%calculate the time vector from tspan and h inputs 
t0 = tspan(1); %initial time
tf = tspan(2); %final time
time = (t0:h:tf)';

n = length(time); %number of elements in the time vector

y = zeros(n,2); %preallocate values for altitude and velocity 

%use a for loop to calculate ODE using the midpoint method 
for i=1:n-1
    %calculate the predicted midpoint value
    ymid = y(i,:) + dydt(time(i),y(i,:))'*h/2;
    %use predicted midpoint value ymid for the slope for next time step
    y(i+1,:) = y(i,:) + dydt(time(i)+h/2,ymid)'*h;
end


altitude = y(:,1); %assign the first column vector as the altitude
velocity = y(:,2); %assign the second column vector as the velocity

end 