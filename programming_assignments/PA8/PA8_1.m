%Example function call used in first tests
yprime = @(t,y) t-y;
t = (0:0.2:2.8)';
[yValues,ea] = PA8_heun(yprime,t,1);

function [y, ea] = PA8_heun(dydt, time, y0, es, varargin)
% [y, ea] = PA8_heun(dydt, time, y0, es, p1, p2, ...)
%
% Heun's Method to solve single ODEs
%
% Inputs:
%   dydt      : An anonymous function that defines the ODE to be solved. 
%               Accomodates optional parameters 
%   time      : An evenly spaced time vector that defines the step size 
%               and time span
%   y0        : The initial condition y0
%   es        : A stopping criterion for the corrector iteration
%   p1,p2,... : additional parameters used by dydt
%
% Outputs:
%   y         : a column vector of y values corresponding to input time vector
%   ea        : a column vector of ea of each time step with corrector iteration

%define the time vector and the last element in the time span
t = time; tf = t(end);
n = length(t); %number of elements in time span
h = t(2)-t(1); %step size, h

%preallocate vectors to y, ea, ypre, and ycor
y=y0*ones(n,1);
ea_v = ones(n,1);
ypre = ones(n,1);
ycor = ones(n,1);


if t(n)<tf  %add one more t value if h increments don't end at tf 
  t(n+1) = tf;
  n = n+1;
end

%If we don't use the corrector iteration for this method
if nargin<4 || isempty(es)

%use for loop to calculate ODE corresponding to each time vector 
for i = 1:n-1
    ea = ones(n,1);
    %calculate the predicted y values
    ypre(i+1) = y(i) + dydt(t(i),y(i),varargin{:})*h;
    %calculate the corrected y values
    ycor(i) = y(i) + ((dydt(t(i),y(i),varargin{:}) + dydt(t(i+1),ypre(i+1),varargin{:}))/(2))*h;
    y(i+1) = ycor(i); %assign the corrected y values to the y vector 
    ea = abs((y(i+1)-ypre(i+1))/(y(i+1))); %calculate absolute value of ea
end

% If we use the predictor-corrector iteration method
else
%use for loop to calculate ODE corresponding to each time vector
for i = 1:n-1
    %calculate the predicted y values 
    ypre(i+1) = y(i) + dydt(t(i),y(i),varargin{:})*h;
    ea = 1; %assign ea
    %do a while loop as long as absolute value of ea is bigger than the stopping criterion
    while abs(ea) > es
        %calculate the corrected y values
        ycor(i+1) = y(i) + ((dydt(t(i),y(i),varargin{:}) + dydt(t(i+1),ypre(i+1),varargin{:}))/(2))*h;
        ea = (ycor(i+1)-ypre(i+1))/(ycor(i+1)); %calculate ea
        ypre(i+1) = ycor(i+1); %assign the corrected y as the new predicted y
    end
    y(i+1) = ycor(i+1); %assign the most recent corrected y values to the y vector 
    ea_v(i+1) = ea; %assign ea to the ea vectors 
end 
    
ea = abs(ea_v); %make all values in ea as absolute values 
ea(1)=[]; %delete the first element of ea since no error associated with y0
end
    
end