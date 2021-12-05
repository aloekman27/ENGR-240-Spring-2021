%Function call for all scalar inputs
R = [333 150 1000 15000 3000 500 2000 100];
V = [80 50 -40 60];
[conditionNumber, currents, dVectors] = circuitSolver(R, V(1), V(2), V(3), V(4));
%Function call when V2 input is a vector
V2vector = 30:0.5:120;
[conditionNumber2, currentsMatrix, dVectorMatrix] = circuitSolver(R, V(1), V2vector, V(3), V(4));
%Note that the condition number of A should not change

function [conditionNumber, currents, dVectors] = circuitSolver(R, V1, V2, V3, V4)
%[conditionNumber, currents, dVectors] = circuitSolver(R, V1, V2, V3, V4)
%
% This function would find us the condition number of the resistor values (R) in the circuit,
% the value of currents, and intermediate vector D using user defined inputs
% 
% Inputs = 
%   R  = Array of resistor values from 1-8. (There must be 8 values because 8 resistors in the circuit) in Ohms
%   V1 = A scalar value for voltage 1. in Volts
%   V2 = A scalar or vector value for voltage 2. in Volts
%   V3 = A scalar value for voltage 3. in Volts
%   V4 = A scalar value for voltage 3. in Volts
%
% Outputs = 
%   conditionNumber = finds the 2-norm condition number of the coefficient matrix A
%   currents = A four column matrix of the current results corresponding to input voltages V1-V4
%   dVectors = A four column intermediate matrix resulting from forward substitution using LU Factorization.


% Work by hand to solve the linear system and build a matrix of the resistor values.
A = [R(6)+R(1)+R(2) -R(1) -R(2) -R(6); R(1) -R(8)-R(3)-R(4)-R(1) R(4) R(8); 
    -R(2) -R(4) R(2)+R(4)+R(5) 0; R(6) R(8) 0 -R(8)-R(6)-R(7)];

%Find the 2-norm condition number of coefficient matrix A
conditionNumber = cond(A);


inverseA = inv(A); %calculate the inverse of matrix A, this will be an intermediate step to find currents
[L,U,P] = lu(A); %do the LU factorization for matrix A. 

%preallocate matrix values for currents and dVector variables. 
currents = ones(length(V2),4);
dVectors = ones(length(V2),4);


%We will use a for loop to test sensitivity of V2 which could be either scalar or vector
for ndx = 1:length(V2)
    % b is equal to the right hand side matrix, which are the voltages
    b = [V1;V2(ndx);V3;V4-V2(ndx)];
    % calculate the Current, I, using matrix multiplication
    I = (inverseA*b)';
    % calculate intermediate vector d using Left Division
    dV = (L\(P*b))';
    
    %Put the values of I and dV into the matrix we preallocated above
    currents(ndx,:) = I;
    dVectors(ndx,:) = dV;
end

end