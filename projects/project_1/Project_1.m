%A script m-file that is used to call the function Mach_Number over several
%conditions
%The conditions will vary based on subsonic and supersonic conditions as
%well as the different values of k

%The two function calls below will use the value k for CO2 in subsonic and
%supersonic conditions (k(CO2)=1.285
tic
[A1,M1] = Mach_Number([0.001 0.002],1.285);
figure(1)
plot(A1,M1,'g') %Subsonic Condition
title('Subsonic Flow for CO2')
xlabel('Cross Sectional Area, A')
ylabel('Mach Number, M')

[A2,M2] = Mach_Number([2 4],1.285);
figure(2)
plot(A2,M2,'g') %Supersonic Condition
title('Supersonic Flow for CO2')
xlabel('Cross Sectional Area, A')
ylabel('Mach Number, M')

%The two function calls below will use the value k for Air in subsonic and
%supersonic conditions (k(Air)=1.400
[A3,M3] = Mach_Number([0.001 0.002],1.4);
figure(3)
plot(A3,M3) %Subsonic Condition
title('Subsonic Flow for Air')
xlabel('Cross Sectional Area, A')
ylabel('Mach Number, M')

[A4,M4] = Mach_Number([2 4],1.4);
figure(4)
plot(A4,M4) %Supersonic Condition
title('Supersonic Flow for Air')
xlabel('Cross Sectional Area, A')
ylabel('Mach Number, M')

%The two function calls below will use the value k for Noble Gasses in subsonic and
%supersonic conditions (k(Noble Gasses)=1.667
[A5,M5] = Mach_Number([0.001 0.002],1.667);
figure(5)
plot(A5,M5,'r') %Subsonic Condition
title('Subsonic Flow for Noble Gasses')
xlabel('Cross Sectional Area, A')
ylabel('Mach Number, M')

[A6,M6] = Mach_Number([2 4],1.667);
figure(6)
plot(A6,M6,'r') %Supersonic Condition
title('Supersonic Flow for Noble Gasses')
xlabel('Cross Sectional Area, A')
ylabel('Mach Number, M')



%plot the graph for subsonic mach values 
figure (7)
plot(A1,M1,A3,M3,A5,M5)
title('Mach Number Values for Supersonic Behavior')
xlabel('Cross Sectional Area, A')
ylabel('Mach Number, M')
legend('k_1 = 1.285 (CO_2)', 'k_2 = 1.400 (Air)', 'k_3 = 1.667 (Noble Gasses)')

%plot the graph for supersonic mach values 
figure (8)
plot(A2,M2,A4,M4,A6,M6)
title('Mach Number Values for Supersonic Behavior')
xlabel('Cross Sectional Area, A')
ylabel('Mach Number, M')
legend('k_1 = 1.285 (CO_2)', 'k_2 = 1.400 (Air)', 'k_3 = 1.667 (Noble Gasses)')
toc

figure (9)
plot(A1,M1,A2,M2)
title('Subsonic and Supersonic Values for CO_2')
xlabel('Cross Sectional Area, A')
ylabel('Mach Number, M')
toc
