%Example function call used in first two tests
[xSol, Jinitial, residualsNorm, iterCount] = PA5_xySolver(1.5, 3, [1;-1], 1e-6);

function [xSol, Jinitial, residualsNorm, iterCount] = PA5_xySolver(A, B, xInit, es)
% [xSol, Jinitial, residualsNorm, iterCount] = PA5_xySolver(A, B, xInit, es)
%
% A function that utilizes Newton-Raphson iteration to solve for x and y given values of
% parameters in a system of equations
%
% Inputs: 
%   A             = A scalar value that will be the parameter of A
%   B             = A scalar value that will be the parameter of B
%   xInit         = A column vector of initial guesses for x_0 and y_0
%   es            = A stopping criterion for the iteration
%
% Outputs: 
%   xSol          = A two-element column vector of the solutions for x and y
%   Jinitial      = The Jacobian matrix evaluated at the initial guesses x_0 and y_0
%   residualsNorm = The Euclidean norm of the residuals vector associated with the solution
%   iterCount     = The scalar value of iterations required for convergence (max iterations = 50)


% Define the anonymous function for the two systems of equations
f = @(x) [exp(A*x(1))+x(2)-0.5*x(2)^2-2.5; B*x(2)-4*cos(3*sqrt(x(1)))];

%set up the vector of initial values x_0 and y_0 to a vector
xi = xInit;

    function [ J ] = J_PA5(x)
        %calculates the Jacobian for the 2x2 from the function
        
        J = ones(2,2); %preallocate a 2x2 matrix
        
        % calculate the partial derivatives of each equation and call it in each position of the matrix
        J(1,1) = A*exp(A*x(1));
        J(1,2) = 1-x(2);
        J(2,1) = (6*sin(3*sqrt(x(1))))/(sqrt(x(1)));
        J(2,2) = B;
    end
    
    % use the function NRsys to help calculate our needed outputs
    [xSol,residualsNorm,ea,iterCount] = NRsys(f,@J_PA5,xi,es,[]);
    
    % set up the values from the function NRsystems and rename them to our desired variable names according to 
    % MATLAB Grader
    xSol
    Jinitial = J_PA5(xi)
    residualsNorm = norm(residualsNorm)
    iterCount

end