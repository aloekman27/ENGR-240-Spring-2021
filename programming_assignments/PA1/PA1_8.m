%input values for v0 and y0
v0 = 25;
y0 = 3.5;

% three different initial angles (in degrees)
theta1 = 30;
theta2 = 45;
theta3 = 60;

%acceleration due to gravity (m/s^2)
g = 9.81;

%three equations that defines d for each initial angle
d1 = ((v0*cosd(theta1))/g) * ((v0*sind(theta1) + sqrt((v0*sind(theta1))^2 + 2*g*y0))); %theta1
d2 = ((v0*cosd(theta2))/g) * ((v0*sind(theta2) + sqrt((v0*sind(theta2))^2 + 2*g*y0))); %theta2
d3 = ((v0*cosd(theta3))/g) * ((v0*sind(theta3) + sqrt((v0*sind(theta3))^2 + 2*g*y0))); %theta3

% x values ranging from 0 to d1/2/3 using data from above with theta1/2/3
x1 = linspace(0,d1,200)'; %theta1
x2 = linspace(0,d2,200)'; %theta2
x3 = linspace(0,d3,200)'; %theta3

% y1/2/3 values that correspond to x1/2/3 due to theta1/2/3
y1 = x1.*tand(theta1) - ((1/2).*(x1.^2.*g)./((v0.*cosd(theta1)).^2)) + y0; %theta1
y2 = x2.*tand(theta2) - ((1/2).*(x2.^2.*g)./((v0.*cosd(theta2)).^2)) + y0; %theta2
y3 = x3.*tand(theta3) - ((1/2).*(x3.^2.*g)./((v0.*cosd(theta3)).^2)) + y0; %theta3

%plot x1/2/3 and y1/2/3 to the same graph
plot (x1,y1,x2,y2,x3,y3)

%annotations for the graph
title('Projectile Trajectories at Different Initial Angles')
xlabel('distance (m)')
ylabel('height (m)')
legend('\theta = 30\circ','\theta = 45\circ','\theta = 60\circ')
grid