%example function call for first three tests
[erfTrapezoid, erfGL3, erfIntegral] = PA7_erfIntegration(1.23)

function [erfTrapezoid, erfGL3, erfIntegral] = PA7_erfIntegration(x)
% [erfTrapezoid, erfGL3, erfIntegral] = PA7_erfIntegration(x)
%
% Use 3 numerical integration approaches to evaluate the error function integral
%       (trapz/cumtrapz, Gauss-Legendre, MATLAB's 'integral')
%
% Input :
%   x            : A column vector of one or more x values at which erf(x) is to be computed
%
% Output :
%   erfTrapezoid : The estimate for erf(x) using composite trapezoid rule
%   erfGL3       : The estimate for erf(x) using three-point Gauss-Legendre quadrature
%   erfIntegral  : The estimate for erf(x) using MATLAB's integral function

%Solves the estimate for erf(x) using the composite trapezoidal rule
%This is if the input x is a scalar
if length(x) == 1
    t = [0, x];
    v = exp(-t.^2);

    I_trapz = trapz(t,v);

%if the x is a vector starting at 0
elseif length(x) ~= 1 && x(1) == 0
    t = x;
    v = exp(-t.^2);
    
    I_trapz = cumtrapz(t,v);

%if x is a vector and does not start at 0
else 
    t = [0; x]; 
    v = exp(-t.^2);
    
    I_trapz = cumtrapz(t,v);
    I_trapz(1)=[]; %deleting the first elemet in the vector because it's not needed in our answer
   
end 

%calculate the value of erf using the trapezoidal rule method
erfTrapezoid = ((2/sqrt(pi)).*I_trapz);


%a new function to define shifted limits of integration
v_new = @(t) (x./2).*((2/sqrt(pi)).*exp(-((x./2)+(x./2).*t).^2));

%calcualte the value of erf using Gauss-Legendre method
erfGL3 = ((5/9.*v_new(-sqrt(3/5))) + (8/9.* v_new(0)) + (5/9.*v_new(sqrt(3/5))));



%define a function of t that will be later used for the integral function 
fun = @(t) exp(-t.^2);
I_Int = zeros(length(x),1);%preallocate vectors for the integral function
for ndx = 1:length(x)
    I_Int(ndx) = integral(fun,0,x(ndx));
end 
%calculate the erf value using MATLAB's integral function
erfIntegral = ((2/sqrt(pi)).*I_Int);
    
    
end 