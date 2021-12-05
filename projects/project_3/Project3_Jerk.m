load project3_accel.dat


t = project3_accel(:,1);
a = project3_accel(:,2);
n = length(t);
h = t(2)-t(1);
dadt_fdf = ones(n,1);


dadt_fdf(1) = (-a(3)+4*a(2)-3*a(1))/(2*h);
dadt_fdf(2) = (-a(4)+4*a(3)-3*a(2))/(2*h);

for k = 3:n-2
    dadt_fdf(k) = (-a(k+2)+8*a(k+1)-8*a(k-1)+a(k-2))/(12*h);
end 

dadt_fdf(n-1) = (3*a(n-1)-4*a(n-2)+a(n-3))/(2*h);
dadt_fdf(n) = (3*a(n)-4*a(n-1)+a(n-2))/(2*h);

xx = min(t):0.01:max(t);
fx_clamped = [0, a', 0];
yy_spline_clamped = spline(t,fx_clamped,xx);
h_new = xx(2)-xx(1);
n_int = length(xx);
dadt_int_fdf = ones(length(xx),1);

dadt_int_fdf(1) = (-yy_spline_clamped(3)+4*yy_spline_clamped(2)-3*yy_spline_clamped(1))/(2*h_new);
dadt_int_fdf(2) = (-yy_spline_clamped(4)+4*yy_spline_clamped(3)-3*yy_spline_clamped(2))/(2*h_new);
                
for p = 3:n_int-2
    dadt_int_fdf(p) = (-yy_spline_clamped(p+2)+8*yy_spline_clamped(p+1)-8*yy_spline_clamped(p-1)+yy_spline_clamped(p-2))/(12*h_new);
end 

dadt_int_fdf(n_int-1) = (3*yy_spline_clamped(n_int-1)-4*yy_spline_clamped(n_int-2)+yy_spline_clamped(n_int-3))/(2*h_new);
dadt_int_fdf(n_int) = (3*yy_spline_clamped(n_int)-4*yy_spline_clamped(n_int-1)+yy_spline_clamped(n_int-2))/(2*h_new);


dadt_grad = gradient(a,h);

dadt_grad_int = gradient(yy_spline_clamped,0.01);


figure(1)
plot(t,a,'rd',xx,yy_spline_clamped,'g-')
title('Acceleration vs. Time Graph')
xlabel('time (s)')
ylabel('acceleration (m/s^{2})')
legend('data','spline interpolation')

figure(2)
plot(xx,yy_spline_clamped,'g-',xx,dadt_grad_int,'r--');
title('Acceleration & Jerk vs. Time Graph using Spline interpolation')
xlabel('time(s)')
ylabel('acceleration(m/s^{2})')
legend('data','jerk(gradient method) m/s^{3}')

figure(3)
plot(t,a,'rd',t,dadt_fdf,'k--',t,dadt_grad,'g--')
title('Acceleration & Jerk vs. Time Graph')
xlabel('time(s)')
ylabel('acceleration(m/s^{2})/jerk(m/s^{3})')
legend('data','finite difference','gradient')

figure(4)
AX = plotyy(t,a,t,dadt_grad);
title('Acceleration & Jerk vs. Time Graph')
xlabel('time (s)')
ylabel(AX(1),'acceleration(m/s^{2})')
ylabel(AX(2),'jerk(m/s^{3})')

figure(5)
AT = plotyy(t,a,t,dadt_fdf);
title('Acceleration & Jerk vs. Time Graph')
xlabel('time (s)')
ylabel(AT(1),'acceleration(m/s^{2})')
ylabel(AT(2),'jerk(m/s^{3})')
legend('acceleration','jerk(finite difference)')
