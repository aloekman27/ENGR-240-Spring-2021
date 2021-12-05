%Function call used for first test
yprime = @(t,y) t-y.^2;
t = (0:0.2:2.8)';
[y, ea] = PA9_odeABM4(yprime,t,1);

function [y, ea] = PA9_odeABM4(dydt,t,y0,es,varargin)
% [y, ea] = PA9_odeABM4(dydt,t,y0,es,varargin)
%
% 4th order Adams-Bashforth-Moulton to solve a differentail equation
%
% Inputs:
%   dydt      : An anonymous function that defines the ODE to be solved. 
%               Accomodates optional parameters 
%   t         : An evenly spaced time vector that defines the step size 
%               and time span
%   y0        : The initial condition y0
%   es        : A stopping criterion for the corrector iteration
%   p1,p2,... : additional parameters used by dydt
%
% Outputs:
%   y         : a column vector of y values corresponding to input time vector
%   ea        : a column vector of ea of each time step with corrector iteration
  
%define the time values for this function
tf = t(end);
n = length(t); %the number of elements in time to be solved
h = t(2)-t(1); %the step size, h

% if necessary, add an additional value of t
% so that range goes from t = ti to tf
if t(n)<tf;
  t(n+1) = tf;
  n = n+1;
end

%Preallocate vectors for the solution
y = y0*ones(n,1); %yValues vector
ypre = ones(n,1); %predicted yValues
ycor = ones(n,1); %corrected yValues
ea_v = ones(n,1); % ea values for each corrector iteration
ea = ones(n,1);
%Use 4th order Runge-Kutta Method for the first three time steps
%       y(i+1) = yi + 1/68(k1+2k2+2k3+k4)*h
for i = 1:3
    k1 = dydt(t(i),y(i),varargin{:})'; %find k1
    ymid2 = y(i) + k1*h/2; %ymid2 used to find k2
    k2 = dydt(t(i)+h/2,ymid2,varargin{:})';  %find k2 
    ymid3 = y(i) + k2*h/2; %ymid3 used to find k3
    k3 = dydt(t(i)+h/2,ymid3,varargin{:})';  %find k3
    yend = y(i) + k3*h;%yend used to find k4
    k4 = dydt(t(i)+h,yend,varargin{:})'; %find k4  
    phi = (k1+2*k2+2*k3+k4)/6; %phi from k1,k2,k3,and k4
    y(i+1) = y(i) + phi*h; %result using 4th orger Runge-Kutta
    ea = (y(i+1)-y(1))/(y(i+1));
end

%if we don't use the corrector iteration
if nargin<4 || isempty(es)
    %calculate using for loop from the 4th time step to the end 
    for i = 4:n-1
        %calculate the ypre values using 4th order Adams-Bashforth 
        ypre(i+1) = y(i) + (h/24)*(55*dydt(t(i),y(i),varargin{:})-59*dydt(t(i-1),y(i-1),varargin{:})+37*dydt(t(i-2),y(i-2),varargin{:})-9*dydt(t(i-3),y(i-3),varargin{:}));
        %calculate the ycor values using 4th order Adams-Moulton
        ycor(i+1) = y(i) + (h/24)*(9*dydt(t(i+1),ypre(i+1),varargin{:})+19*dydt(t(i),y(i),varargin{:})-5*dydt(t(i-1),y(i-1),varargin{:})+dydt(t(i-2),y(i-2),varargin{:}));
        %set the corrected y to the vector y
        y(i+1) = ycor(i+1);
    end 

%if we use the corrector iteration
else 
    %calculate using for loop from the 4th time step to the end 
    for i = 4:n-1
        %calculate the ypre values using 4th order Adams-Bashforth 
        ypre(i+1) = y(i) + (h/24)*(55*dydt(t(i),y(i),varargin{:})-59*dydt(t(i-1),y(i-1),varargin{:})+37*dydt(t(i-2),y(i-2),varargin{:})-9*dydt(t(i-3),y(i-3),varargin{:}));
        ea = 1; %preallocate ea 
        
        %do a while loop as long as absolute value of ea is bigger than the stopping criterion
        while abs(ea) > es
            %calculate the ycor values using 4th order Adams-Moulton
            ycor(i+1) = y(i) + (h/24)*(9*dydt(t(i+1),ypre(i+1),varargin{:})+19*dydt(t(i),y(i),varargin{:})-5*dydt(t(i-1),y(i-1),varargin{:})+dydt(t(i-2),y(i-2),varargin{:}));
            %calculate the ea of the iteration
            ea = (ycor(i+1)-ypre(i+1))/(ycor(i+1));
            %update ypre values
            ypre(i+1) = ycor(i+1);
        end 
        %set the corrected y to the vector y
        y(i+1) = ycor(i+1);
        ea_v(i+1) = ea; %assign ea to the vector
    end
    %absolute value of the error
    ea = abs(ea_v);
    ea(1:4)=[]; %delete elements 1:4 of ea because it is not included in the iteration
    
end
end