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


tic
dadt_fit_grad = gradient(amodel(tm,A,B,C,D),0.01);
toc

figure(2)
AP=plotyy(tm,amodel(tm,A,B,C,D),tm,dadt_fit_grad);
title('Acceleration & Jerk vs. Time Curve Fit Differentiation')
xlabel('time (s)')
ylabel(AP(1),'acceleration (m/s^{2})')
ylabel(AP(2),'jerk (m/{s^3})')
legend('curve fit','gradient differentiation')








