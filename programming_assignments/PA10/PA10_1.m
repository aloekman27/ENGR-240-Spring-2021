%function call used in first two tests
[hLinear, hNonlinear] = PA10_aquifer((0:50:1000)',[12 5], 1.2, 1.5e-4);

function [hLinear, hNonlinear] = PA10_aquifer(xspan, h_BC, K, N)
% [hLinear, hNonlinear] = PA10_aquifer(xspan, h_BC, K, N)
%
% Finding the height of an unconfined groundwater aquifer using non-linear and
% linear differential equations solved using Boundary Values
%
% Inputs :
% xspan : A column vector of x values at which h(x) is to be calculated
% h_BC : A two-element vector of the boundary conditions
% K : A scalar input for the hydraulic permeability
% N : A scalar input for the infiltration rate
%
% Outputs :
% hLinear : h(x) values corresponding to x with simplified linear ODE
% hNonlinear : h(x) values corresponding to x with nonLinear ODE model

%Start of Linear ODE
%define the boundary conditions (meters)
h0 = h_BC(1);
hL = h_BC(2);

%average water table height used for linear ODE (meters)
hbar = (h0+hL)/2;

%define an anonymous function for the linear model
yp = @(x,y) [y(2);-N/(hbar*K)];

guesses = [-100 100]; %initial guess
[x1,y1] = ode45(yp,xspan,[h0 guesses(1)]); %use ode45 with first guess
shot1 = y1(end,1); %save the last point (boundary condition) 
[x2,y2] = ode45(yp,xspan,[h0 guesses(2)]); %use ode45 with second guess
shot2 = y2(end,1); 
shots = [shot1 shot2]; %input the last points to a vector 
newIC = interp1(shots,guesses,hL,'linear','extrap');%interpolate the last point with the height
[x3,y3] = ode45(yp,xspan,[h0 newIC]); %use the newIC as the initial value 
hLinear = y3(:,1); %the first column outputs all h values for the corresponding x values


figure(1)
plot(x3,hLinear)
title('Solution of Linear Differential Equation')
xlabel('x (m)')
ylabel('h (m)')

%End of Linear ODE
%Start of NonLinear ODE

%define a function to calculate the height from the nonlinear model
    function height = height_1(xspan,h_BC,K,N)
    % height = height_1(xspan,h_BC,K,N)
    % 
    % a function for the nonlinear model to calculate the height from the 
    % corresponding x values
    % 
    % Inputs :
    %   xspan : A column vector of x values at which h(x) is to be calculated
    %   h_BC  : A two-element vector of the boundary conditions
    %   K     : A scalar input for the hydraulic permeability
    %   N     : A scalar input for the infiltration rate
    %
    % Outputs :
    %   height : height of the water table from the corresponding x values using
    %            nonlinear ODE model

    L = xspan(end); %the last x value on the interval xspan
    h0 = h_BC(1); %boundary condition on x = xspan(1)
    hL = h_BC(2);%boundary condition on x = xspan(end)

    %an anonymous function to define the nonlinear ODE model
    ypnonl = @(x,y) [y(2);-N/(y(1)*K)-(y(2)^2/(y(1)))];

        %a function that will be used to find the roots of nonlinear model
        function r = bar_res(IC_guess)
        % r = bar_res(IC_guess)
        %
        % Find the roots, r, of a nonlinear model
        % 
        % Inputs :
        %   IC_guess = guess for the inital value for the ODE
        %
        % Outputs : 
        %   r = roots of the nonlinear model

        %use ode45 to solve the nonlinear model
        [~,y] = ode45(ypnonl,xspan,[h0 IC_guess]);
        %find the roots from the ODE solver
        r = y(end,1)-hL;

        end
        
     % find the right initial value using the bar_res function
    ICguess_on_target = fzero(@(x) bar_res(x),-0.0001);
    %solve the ODE using ode45
    [~,y] = ode45(ypnonl,xspan,[h0 ICguess_on_target]);
    % the first column is the height corresponding to x values
    height = y(:,1);
   
    end
    %use the height_1 function to solve for the height from the nonlinear model    
hNonlinear = height_1(xspan,[h0 hL],K,N);

figure(2)
plot(xspan,hNonlinear)
title('Solution of non-Linear Differential Equation')
xlabel('x (m)')
ylabel('h (m)')
end