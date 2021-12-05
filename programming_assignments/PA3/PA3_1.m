%This is the function call used for the first three tests
f = @(x) sin(2*x) - tan(x/4);
[xr,res,ea,iter] = PA3_secant(f, [1 1.1], 1e-3, 50);

function [xroot,residual,ea,iterCount] = PA3_secant(fun,xi,es,max_it,varargin)
%[xroot,residual,ea,iterCount] = PA3_secant(fun,xi,es,max_it,p1,p2,....)
%The Secant Method for the solution of a roots problem

%Inputs
%fun : the function defined by the user that will be solved using the numerical method to find the roots
%xi : initial guesses, this should be a 1,2 vector where the 1st vector is the 1st guess and 2nd vector is the 2nd guess
%es : optional stopping criterion, default is es = 1E-6
%max_it :optional maximum iterations, default is max_it = 30
%pi,p2,... : optional parameter inputs for function evaluation

%Outputs
%xroot : root estimate 
%residual : the function's value at x = xroot (root estimate)
%ea : approximate relative error in root estimate
%iterCount : number of iterations required for convergence


%copy/paste the commands for your function here.  You can change the names, 
%but not the order, of the input and/or output variables in the function 
%command to match your code.  Do not change the name of the function, 
%however, as this exact name (i.e. "PA3_secant") is required for 
%the test suite to run.

%set up each element of the input guesses into two different variables
xroot1 = xi(1); xroot2 = xi(2);
%Determine errors and default values if user input is insufficient for the code to run
if nargin < 2, error('Two few inputs, read help comments'), end %not enough input
if nargin < 3||isempty(es), es=1e-5; end %stopping criterion default
if nargin < 4||isempty(max_it), max_it=30; end %maximum iterations default


%Set up a for loop to determine the xroot of the function until it converges
for iterCount = 1:max_it
    %The formula that is used to solve for xroot using the secant method
    xroot = xroot2-(fun(xroot2,varargin{:})*(xroot1-xroot2))/(fun(xroot1,varargin{:})-fun(xroot2,varargin{:})); 
    %the second guess is now the first guess in the following loop
    xroot1 = xroot2;
    %the xroot we found is now the second guess
    xroot2 = xroot;
    
    %calculate approximate error.
    if xroot ~= 0, ea = (xroot-xroot1)/xroot; end
    
    if abs(ea) <= es, break, end 
end 
%calculate residual value of the function
residual = fun(xroot,varargin{:});

end