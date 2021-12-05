% enter known values
R = 8.31447; % ideal gas constant (J/mol-K)
MO2 = 15.999; % Molar Mass (g/mol)
MCO2 = 44.01; % Molar Mass (g/mol)
V = 0.15; % Volume (m^3)

m = 100; %mass of gas (g)

%create 100 equally spaced values from the question above
P = linspace(500*10^3,100*10^3)';

%rewrite the equation so we can find the dependent variable T
TO2 = (P.*V.*MO2)./(m.*R) %for O2 gas

TCO2 = (P.*V.*MCO2)./(m.*R) %for CO2 gas
plot  (P, TO2)