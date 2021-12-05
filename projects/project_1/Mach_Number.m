function [A,M_value] = Mach_Number(xi,k)

%function [M] = Mach_Number(,xi,es,max_it)
%
%Using the Secant Method to find the value of the Mach Number (M)
%
%Inputs
%   xi:     initial guesses
%   es:     optional stopping criterion, default is es = 1E-5
%   max_it: optional maximum iterations, default is max_it = 30
%
%Outputs
%   M: Mach Number estimate

%Defining the variables
A = (1.0:0.1:10)';
es = 1e-6;
max_it = 50;

%Preallocating the value of M
M_value = zeros(length(A),1);

%Creating a for loop to run through every element of A
for i = 1:numel(A)
    fun = @(M) 1/M * sqrt((2+(k-1)*M^2)/(k+1))^((k+1)/(k-1))-A(i); %Function to solve for M
    xroot1 = xi(1); %First initial guess
    xroot2 = xi(2); %Second initial guess
    %Creating another for loop to solve for the value of M
    for j = 1:max_it
        xroot = xroot2 - fun(xroot2)*(xroot1-xroot2)/(fun(xroot1)-fun(xroot2));
        if xroot ~= 0
            ea = (xroot-xroot2)/xroot; 
        end
        if abs(ea) <= es
            break 
        end
        xroot1 = xroot2;
        xroot2 = xroot;
    end
    M_value(i) = xroot;
end
end