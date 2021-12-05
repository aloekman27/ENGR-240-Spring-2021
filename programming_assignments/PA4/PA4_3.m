%Example code to call the function
diagonals = [1, -2, 1];
b = [125; 200; 100; 150];
[x, diagDominance, normResiduals] = triDiagSolver(diagonals,b);

function [x, diagDominance, normResiduals] = triDiagSolver(diagValues, RHSvector)
%[x, diagDominance, normResiduals] = triDiagSolver(diagValues, RHSvector)
%
% This function m-file accepts user defined inputs that solves a tridiagonal matrix using
% Gauss-Seidel method. 
%
% Inputs:
%   diagValues = a row vector of three numbers that gives the values for the three diagonals respectively
%   RHSvector = A column vector to define the RHS ("b" vector of the linear system)
%
% Outputs:
%   x = Column vector of the solutions to the system of equations
%   diagDominance = A scalar equal to 1 if coefficient matrix is diagonally dominant and 0 if it is not
%   normResidual  = A scalar that is equal to the Euclidean norm of the residuals vector associated with the solution

%the built-in function spadiag to create the desired coefficient matrix
n = length(RHSvector);
e = ones(n,1);
A = spdiags([diagValues(1)*e diagValues(2)*e diagValues(3)*e],-1:1,n,n);
A = full(A)

%preallocate values for the stopping criterion and maximum iterations
es = 1e-5;
max_it = 50;

%preallocate vectors for the initial x guess and so on.
xi = zeros(length(RHSvector),1);

%set intermediate vector d that will be used in Gauss-Seidel method
C = A;
d = RHSvector./diag(A);


%Start of Gauss-Seidel method. 
%Note: I could not call the GS.m function so I copied from the week 4 mfiles, pasted, and modified the code
for i = 1:n  
  C(i,1:n) = A(i,1:n)/A(i,i);
  C(i,i) = 0;
end

%set values for the approximate error and value of x
ea = ones(n,1);
x = xi;

for iter = 1:max_it
    x_old = x;
    for i = 1:n
        x(i) = d(i)-C(i,:)*x;
    end
    if all(x ~= 0), ea = (x - x_old)./x; end
    if max(abs(ea)) <= es, break, end 
end

% Chek diagonal dominance of the coefficient matrix, A
diagDominance = abs(diagValues(2))>(abs(diagValues(1)+abs(diagValues(2))));

%Calculate the residual vector that will be used to calculate the norm of the matrix
residual = A*x-RHSvector;
normResiduals = norm(residual);

end