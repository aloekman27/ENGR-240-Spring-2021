%Example function call used in first test
[V_mp, V_mp_residual] = PA3_solarPower(298, 0.5, 0.5, 1e-6);

function [V_mp, V_mp_residual] = PA3_solarPower(T, V_OC, V_mp_guess, es)
%[V_mp, V_mp_residual] = PA3_solarPower(T, V_OC, V_mp_guess, es)
%This function uses fzero method from user inputs to determine the maximum voltage the solar cell produces

%Inputs:
%T = Temperature (in Kelvin). User is able to write 1 temperature value or a vector of temperature values
%V_OC = The open circuit voltage (Volts)
%V_mp_guess = The user's initial guess of the voltages of the maximum output (volts)
%es = the stopping criterion for the calculation. 

%Output:
%V_mp = the voltage of the maximum output true value
%V_mp_residual = residual values associated with the numerical solution for each V_mp value

%copy/paste the commands from your function here.  You can change the names, 
%but not the order, of the input and/or output variables in the function 
%command to match your code.  Do not change the name of the function, however, 
%as this exact name (i.e. "PA3_solarPower") is required for the test suite to run.

q = 1.6022e-19; % charge on an electron (Coulombs)
kB = 1.3806e-23; % Boltzmann constant (Joules /Kelvin)

% using the function input, es, as the stopping criterion in the fzero method
options = optimset('TolX',es); 

for ndx = 1 : numel(T)
    %define the function that we want to determine
    fun = @(Vmp) exp(q.*Vmp/(kB*T(ndx)))*(1+(q.*Vmp/(kB*T(ndx))))-exp(q*V_OC/(kB*T(ndx)));
    %define anonymous function with new value of T for bisection method
    [Vmp,residual] = fzero(@(x) fun(x),V_mp_guess, options);
    
    %use the vecotrs below to put the values of V_mp for values of T that are vectors
    V_mp(ndx,1) = Vmp
    V_mp_residual(ndx,1) = residual


end
    
    
    

end