[x, y, xyPeak] = projectileTrajectory(20, 35, 5 )

function [x, y, xyPeak] = projectileTrajectory( v0, theta, y0 )
%Enter the commands for your function in the space below.  Use the variable names specified in 
%the function definition command in line 1.  Be sure to assign values to each of the output variables.

%constant g
g = 9.81; % m/s^2

%rewrite equation to find d in MATLAB code
d = ((v0*cosd(theta))/g) * ((v0*sind(theta) + sqrt((v0*sind(theta))^2 + 2*g*y0))); % d = 44.4686m

%create the array of 200 x values rangin from 0 to d
x = linspace(0,d,200)';

%Find each y values for every x ranging from 0 to d
y = x.*tand(theta) - ((1/2).*(x.^2.*g)./((v0.*cosd(theta)).^2)) + y0;

%find the maximum point of y. p is the assigned variable for the x value when y is max
[maxy,p] = max(y);

%x(p) is the value of x when y is max (max(y))
maxx = x(p);

% the array of coordinates where the trajectory is at its peak
xyPeak = [maxx, maxy];
 
 
plot(x,y)

end