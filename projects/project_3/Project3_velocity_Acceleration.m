load project3_accel.dat % loads the data for the project


t = project3_accel(:,1);
a = project3_accel(:,2);

tic
I_trapz = cumtrapz(t,a);

v0 = -1.6;

vt = v0 + I_trapz;
toc

tic
I_trapz_2 = cumtrapz(t,vt);

s0 = 0;

st = s0 + I_trapz_2;
toc

tic
xx = min(t):0.01:max(t);
fx_clamped = [0, a', 0];
yy_spline_clamped = spline(t,fx_clamped,xx);


I_spline_v = cumtrapz(xx,yy_spline_clamped);

vt_spline = v0 + I_spline_v;
toc

tic
I_spline_s = cumtrapz(xx,vt_spline);

st_spline = s0 + I_spline_s;
toc

figure(1)
AR = plotyy(t,a,t,vt);
title('Acceleration & Velocity vs. Time Integral')
xlabel('time (s)')
ylabel(AR(1),'acceleration (m/s^{2})')
ylabel(AR(2),'velocity (m/s)')

figure(2)
BB = plotyy(t,vt,t,st);
title('Velocity & Position vs. Time Integral')
xlabel('time (s)')
ylabel(BB(1),'velocity (m/s)')
ylabel(BB(2),'position (m)')

figure(3)
plot(t,a,'k-',t,vt,'g-',t,st,'r-')
title('Acceleration, Velocity, and Position vs. Time')
xlabel('time (s)')
legend('acceleration (m/s^{2})','velocity (m/s)','position (m)')

figure(4)
plot(t,a,'rd',xx,yy_spline_clamped,'g-')
title('Acceleration vs. Time Graph')
xlabel('time (s)')
ylabel('acceleration (m/s^{2})')
legend('data','spline interpolation')

figure(5)
CC = plotyy(xx,yy_spline_clamped,xx,vt_spline);
title('Velocity using Interpolation of Acceleration')
xlabel('time (s)')
ylabel(CC(1),'acceleration(m/s^{2})')
ylabel(CC(2),'velocity(m/s)')

figure(6)
DD = plotyy(xx,vt_spline,xx,st_spline);
title('Position using Interpolation of Acceleration')
xlabel('time (s)')
ylabel(DD(1),'velocity(m/s)')
ylabel(DD(2),'position(m)')






