%Calling command for first tests
[coeffMatrix, RHS, r, Temp] = PA10_heat([100 20],[1 3] ,0.2,50.0);

function [coeffMatrix, RHS, r, Temp] = PA10_heat(TBoundaries, radiusInnerOuter, deltaR, Source)
% [coeffMatrix, RHS, r, Temp] = PA10_heat(TBoundaries, radiusInnerOuter, deltaR, Source)
%
% Using Finite Difference Method to solve a Radial Conduction Heat Transfer Differential equation
% with second-order centered finite difference formulas
%
% Input :
%   TBoundaries         : A two-element vector defining the boundary conditions for the solution
%   radiusInnerOuter    : A two-element vector specifying the inner and outer radius of the hollow
%                         cylinder
%   deltaR              : A scalar input for the node spacing
%   Source              : A scalar input for the heat source term
%
% Output :
%   coeffMatrix         : The coefficient matrix for the linear system of equations
%   RHS                 : The right hand side vector for the linear system of equations
%   r                   : A column vector of node locations
%   Temp                : A column vector of the corresponding temperatures at each node in the solution

S = Source; %heat source (Celcius/mm^2)

%define the boundary conditions (Celcius)
y1 = TBoundaries(1); %y at r1 
yn = TBoundaries(end); %y at rn

%define the inner and outer radius used in boundary conditions (mm)
r1 = radiusInnerOuter(1); 
rn = radiusInnerOuter(end);

% define the interval with corresponding nodes
r = (r1:deltaR:rn)';
dr = deltaR;

%calculate the number of nodes 
nodes = numel(r);
n_mat = nodes-2; %internal nodes

%preallocate values for the diagonal in the coefficient matrix
diag1 = ones(n_mat,1);
diagmid = -4*ones(n_mat,1); 
diag2 = ones(n_mat,1);

%use for loop to determine the elements in the first diagonal
for ndx = 2:nodes-1
    diag1(ndx-1) = 2-(dr/r(ndx+1));
end 

%use for loop to determine the elements in the second diagonal
for i=2:nodes-1
    diag2(i-1) = 2+(dr/r(i-1)); 
end

%set the diagonal values into a matrix
diagVal = [diag1 diagmid diag2];
%use spdiags to form the coefficient matrix
coeffMatrix = spdiags(diagVal,-1:1,n_mat,n_mat);

%Calculate the right hand side vector
RHS1 = -2.*S.*(dr^2)-y1.*(2-(dr./r(2))); %RHS 1st element
RHS_int = (-2.*S.*(dr.^2)).*ones(nodes-4,1);%RHS 1st to 2nd to last element
RHSn = -2.*S.*(dr^2)-yn.*(2+(dr./r(nodes-1)));%RHS last element 

%form the RHS vector from the data above
RHS = [RHS1;RHS_int;RHSn]; 

%calculate the yvalues for the internal nodes
y_int = coeffMatrix\RHS;

%calculate the result by combining the boundary conditions
Temp = [y1; y_int; yn];

end