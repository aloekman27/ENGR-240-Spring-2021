%enter parameter values
Ea = 1.11 * 10^5; % Activation energy for the reaction (Joules/mole)
A = 2.126 * 10 ^9;% pre-exponential factor for the reaction 
R = 8.31447; % ideal gas constant (Joules/mole-Kelvin)

%input operator for temperature
tempKelvin = (250:2:400)'; % The column vectors of T (kelvin)

%input equation that will be calculated
reactionRate = A.*exp(-Ea./(R.*tempKelvin)); % k

%Plot the graph using semilog, then annotate the graph
semilogy(tempKelvin,reactionRate) % use semilog for y-axis since it involves exponentials
% Annotations for graph
title('Arrhenius Equation')  
xlabel('Temperature (K)') 
ylabel('Reaction Rate') 
grid
legend('Reaction Rate of Reaction')