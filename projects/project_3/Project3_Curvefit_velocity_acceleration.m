load project3_accel.dat

t = project3_accel(:,1);
a = project3_accel(:,2);


tic
amodel = @(t,A,B,C,D) A.*exp(B.*t).*cos((C.*t)+D);

ares = @(vars,ti,ai) sum((ai-amodel(ti,vars(1),vars(2),vars(3),vars(4))).^2);

options = optimset('TolX',1e-6,'TolFun',1e-6);
vars = fminsearch(ares,[1, 1, 1, 1],options,t,a);

A = vars(1);
B = vars(2);
C = vars(3);
D = vars(4);

St_nonl = sum((a-mean(a)).^2);
Sr_nonl = sum((a-amodel(t,A,B,C,D)).^2);
r_squared_nonl = 1-Sr_nonl/St_nonl
StandardError_nonl = sqrt(Sr_nonl/(length(t)-4))

tm = min(t):0.01:max(t);
toc

figure(1)
plot(t,a,'kd',tm,amodel(tm,A,B,C,D),'r--')
title('Acceleration vs. Time Curve Fit Graph')
xlabel('time (s)')
ylabel('acceleration(m/s^{2})')
legend('data','fminsearch curve fit')


v0 = -1.6;
s0 = 0;
tic
I_fit_v = cumtrapz(tm,amodel(tm,A,B,C,D));
v_fit = v0 + I_fit_v;
toc

figure(2)
EE=plotyy(tm,amodel(tm,A,B,C,D),tm,v_fit);
title('Acceleration & Velocity vs. Time Curve Fit')
xlabel('time (s)')
ylabel(EE(1),'acceleration (m/s^{2})')
ylabel(EE(2),'velocity (m/s)')

tic
I_fit_s = cumtrapz(tm,v_fit);
s_fit = s0 + I_fit_s;

figure(3)
FF=plotyy(tm,v_fit,tm,s_fit);
title('Velocity & Positionvs. Time Curve Fit')
xlabel('time (s)')
ylabel(FF(1),'velocity (m/s)')
ylabel(FF(2),'position (m)')






