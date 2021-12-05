%%              ENGR 240: Applied Numerical Methods
%                           Project 3
% This project will cover the methods of integration and differentiation to
% find the position, velocity, and jerk values using from a damped
% oscillation function and acceleration data, and plot all of them on
% different graphs to compare its values.



load project3_accel.dat %load the data for the project


t = project3_accel(:,1); % time column vector (s)
a = project3_accel(:,2); % acceleration column vector (m/s^2)
n = length(t); %length of t (or a)
h = t(2)-t(1); % find the step of the interval
%% Calculate Jerk using Finite-Difference Formula
% This section will calculate the jerk of the damped oscillation function
% using finite-difference differentiation method from the original data

tic
%Preallocate values for each derivative of the data
dadt_fdf = ones(n,1);

%use forward finite-difference formula for the first two derivatives
dadt_fdf(1) = (-a(3)+4*a(2)-3*a(1))/(2*h);
dadt_fdf(2) = (-a(4)+4*a(3)-3*a(2))/(2*h);

%use centered finite-difference formula for the 3rd and 2nd to last
%derivatives
for k = 3:n-2
    dadt_fdf(k) = (-a(k+2)+8*a(k+1)-8*a(k-1)+a(k-2))/(12*h);
end 

%use backward finite-difference formula for the 2 last derivatives
dadt_fdf(n-1) = (3*a(n-1)-4*a(n-2)+a(n-3))/(2*h);
dadt_fdf(n) = (3*a(n)-4*a(n-1)+a(n-2))/(2*h);
toc

%% SPLINE INTERPOLATION
% This section will compute the jerk value of the damped oscillation using
% a finite-difference differentiation method from an interpolated data
% calculated using a spline function with clamped ends

tic
xx = min(t):0.01:max(t);%calculate the new interval for the time
fx_clamped = [0, a', 0];%clamped end condition
yy_spline_clamped = spline(t,fx_clamped,xx);%spline
h_new = xx(2)-xx(1); %step for the interval
n_int = length(xx); %calculate the length of the time elements 
dadt_int_fdf = ones(length(xx),1); %preallocate a vector for the derivative

%use forward finite-difference formula for the first two derivatives
dadt_int_fdf(1) = (-yy_spline_clamped(3)+4*yy_spline_clamped(2)-3*yy_spline_clamped(1))/(2*h_new);
dadt_int_fdf(2) = (-yy_spline_clamped(4)+4*yy_spline_clamped(3)-3*yy_spline_clamped(2))/(2*h_new);
        
%use centered finite-difference formula for the 3rd and 2nd to last
%derivatives
for p = 3:n_int-2
    dadt_int_fdf(p) = (-yy_spline_clamped(p+2)+8*yy_spline_clamped(p+1)-8*yy_spline_clamped(p-1)+yy_spline_clamped(p-2))/(12*h_new);
end 

%use backward finite-difference formula for the 2 last derivatives
dadt_int_fdf(n_int-1) = (3*yy_spline_clamped(n_int-1)-4*yy_spline_clamped(n_int-2)+yy_spline_clamped(n_int-3))/(2*h_new);
dadt_int_fdf(n_int) = (3*yy_spline_clamped(n_int)-4*yy_spline_clamped(n_int-1)+yy_spline_clamped(n_int-2))/(2*h_new);
toc

tic
%Use the gradient derivative function to calculate the derivative of the
%data
dadt_grad = gradient(a,h);
toc

tic
%Use the gradient derivative function to calculate the derivative of the
%interpolated data
dadt_grad_int = gradient(yy_spline_clamped,0.01);
toc


%% CURVE FIT OF THE DATA
% This section will compute the curved-fitted acceleration data using the
% built-in function "fminsearch" along with the coefficient of
% determination and standard error

tic
%define an anonymous function for the fminsearch
amodel = @(t,A,B,C,D) A.*exp(B.*t).*cos((C.*t)+D);

%calcualte the squared sum of the residuals for fminsearch
ares = @(vars,ti,ai) sum((ai-amodel(ti,vars(1),vars(2),vars(3),vars(4))).^2);

options = optimset('TolX',1e-6,'TolFun',1e-6);
%calculate each variable A,B,C,D for the expected curve fit
vars = fminsearch(ares,[1, 1, 1, 1],options,t,a);

A = vars(1);%variable A
B = vars(2);%variable B
C = vars(3);%variable C
D = vars(4);%variable D

%Calculate the error of the fit
St_nonl = sum((a-mean(a)).^2);
Sr_nonl = sum((a-amodel(t,A,B,C,D)).^2);
r_squared_nonl = 1-Sr_nonl/St_nonl
StandardError_nonl = sqrt(Sr_nonl/(length(t)-4))

tm = min(t):0.01:max(t);
toc

%% DERIVATIVE OF THE CURVE FIT DATA
% This section will compute the jerk values of the damped oscillation using
% a built-in differentiation function "gradient" from the curved-fitted
% data from the previous section

tic
%calculate the derivative of the data using gradient method
dadt_fit_grad = gradient(amodel(tm,A,B,C,D),0.01);
toc

%% INTEGRATION OF THE DATA
% This section will compute the velocity and position values using a
% built-in integration method "cumtrapz" from the original acceleration
% data


tic
% calculate the integral using composite trapezoidal rule to find velocity
I_trapz = cumtrapz(t,a);

v0 = -1.6;%value of v0 (m/s)

vt = v0 + I_trapz;%calculate v(t) from the integration and v0
toc

tic
% calculate the integral using composite trapezoidal rule to find position
I_trapz_2 = cumtrapz(t,vt);

s0 = 0;%value of s0 (m)

st = s0 + I_trapz_2;%calculate s(t) from the integration and s0
toc

%% INTEGRATION OF THE INTERPOLATED DATA
% This section will compute the position and velocity values using a
% built-in integration method "cumtrapz" from the interpolated data

% calculate the integral of the interpolated data using composite
% trapezoidal rule to find velocity
I_spline_v = cumtrapz(xx,yy_spline_clamped);

vt_spline = v0 + I_spline_v;%calculate v(t) of the interpolated data from the 
                            %integration and v0 
toc

tic
% calculate the integral of the interpolated data using composite
% trapezoidal rule to find position
I_spline_s = cumtrapz(xx,vt_spline);

st_spline = s0 + I_spline_s;%calculate s(t) of the interpolated data from the 
                            %integration and s0 
toc

%% INTEGRATION OF THE CURVE FIT
% This section will compute the position and velocity values using a
% built-in integration method "cumtrapz" from the interpolated data

tic
%Calculate the integral of the curve fitted data to find velocity
I_fit_v = cumtrapz(tm,amodel(tm,A,B,C,D));
v_fit = v0 + I_fit_v; %the velocity data of the curve fit 
toc

tic
%Calculate the integral of the curve fitted data to find position
I_fit_s = cumtrapz(tm,v_fit);
s_fit = s0 + I_fit_s;%the position data of the curve fit 
toc

%% PLOT THE GRAPH OF THE RESULTS
% This section will plot of the data that are calculated in 5 separate
% graphs:

%plot the graph of acceleartion and its derivative, jerk
figure(1)
plot(t,a,xx,yy_spline_clamped,tm,amodel(tm,A,B,C,D),t,dadt_fdf,t,dadt_grad)
title('Acceleration & Jerk vs Time Graphs')
xlabel('time (s)')
ylabel('Acceleration (m/s^{2}) / Jerk (m/s^{3})')
legend('data','interpolated data', 'curve fit','jerk(finite difference)','jerk(gradient)')

%plot the graph of acceleartion and its integral, velocity and position
figure(2)
plot(t,a,xx,yy_spline_clamped,tm,amodel(tm,A,B,C,D),t,vt,t,st,xx,vt_spline,xx,st_spline)
title('Acceleration, Velocity, Position vs Time Graphs')
xlabel('time (s)')
ylabel('Acceleration (m/s^{2}) / Velocity(m/s) / Position(m)')
legend('data','interpolated data','curve fit','velocity','"position','velocity(spline)','position(spline)')

%graph of acceleration, velocity, position, and jerk
figure(3)
plot(t,a,t,vt,t,st,t,dadt_fdf)
title('Acceleration, Velocity, Position, and Jerk Graph')
xlabel('time(s)')
legend('acceleration','velocity','"position','jerk')

%graph of acceleration, velocity, position, and jerk from interpolated data
figure(4)
plot(xx,yy_spline_clamped,xx,vt_spline,xx,st_spline,xx,dadt_int_fdf)
title('Acceleration, Velocity, Position, and Jerk Graph for Interpolated Data')
xlabel('time(s)')
legend('acceleration','velocity','"position','jerk')

%graph of acceleration, velocity, position, and jerk from curve fit
figure(5)
plot(tm,amodel(tm,A,B,C,D),tm,v_fit,tm,s_fit,tm,dadt_fit_grad)
title('Acceleration, Velocity, Position, and Jerk Graph for Curve Fit')
xlabel('time(s)')
legend('acceleration','velocity','"position','jerk')


