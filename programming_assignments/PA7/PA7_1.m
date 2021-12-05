load beam_data.dat; %example data file is available on Canvas
%example function call
[M_Oh4FD, M_gradient, M_gradSpline] = PA7_beam(beam_data);

function [M_Oh4FD, M_gradient, M_gradSpline] = PA7_beam(beamDeflectionData)
% [M_Oh4FD, M_gradient, M_gradSpline] = PA7_beam(beamDeflectionData)
% 
% Calculate the second derivative from the data using various methods that will 
% be used to calculate internal bending moment in a beam with the Equation:
%                       
%                   M = E * I * (second derivative of data)
%
% Input :
%   beamDeflectionData : a two column matrix with evenly spaced data for axial 
%                        position, x, in the first column and the corresponding 
%                        lateral deflection measurements, v, in the second column
%
% Output :
%   M_Oh4FD      : The bending moment calculated using 4th order centered finite difference formulas
%   M_gradient   : The bending moment calculated using the gradient function
%   M_gradSpline : The bending moment calculated using the gradient function on the clamped spline 
%                  Interpolation

%Define variables for known data
E = 29*10^6; %psi
I = 156; %in^4
x = beamDeflectionData(:,1); %axial position
v = beamDeflectionData(:,2); %lateral deflection
h = x(2)-x(1); %calculate the step (spacing between x data)

%preallocate a column vector for the second derivative using the first method 
d2vdx2 = ones(length(x),1);

%for loop to calculate the first 2 data using 2nd order forward difference formula
for ndx = 1:2
    d2vdx2(ndx) = (-v(ndx+3)+(4*v(ndx+2))-(5*v(ndx+1))+(2*v(ndx)))/(h^2);
end 

%for loop to calculate the 3rd to 3nd to last data using 4th order center difference formula
for ndx = 3:(length(v)-2)
    d2vdx2(ndx) = (-v(ndx+2)+(16*v(ndx+1))-(30*v(ndx))+(16*v(ndx-1))-v(ndx-2))/(12*(h^2));
end 

%for loop to calculate the last 2 data using 2nd order backward difference formula
for ndx = length(v)-1:length(v)
    d2vdx2(ndx) = ((2*v(ndx))-(5*v(ndx-1))+(4*v(ndx-2))-v(ndx-3))/(h^2);
    
end 

%Calculate the bending moment using the formula given and the second 
%derivative from the method above
M_Oh4FD = E*I.*d2vdx2;

%Use gradient function 2 times to calculate the second derivative of the data
gradient_1 = gradient(v,12);
gradient_2 = gradient(gradient_1,12);

%Calculate the bending moment using the formula given and the second 
%derivative from the method above
M_gradient = E*I.*gradient_2;

%define an interval and a clamped spline for the new data
xx = min(x):0.5:max(x);
fx_clamped = [0, v', 0];
yy_spline_clamped = spline(x,fx_clamped,xx);

%Use gradient function 2 times to calculate the second derivative of the data
spline_gradient_1 = gradient(yy_spline_clamped,0.5);
spline_gradient_2 = gradient(spline_gradient_1,0.5);

%Calculate the bending moment using the formula given and the second 
%derivative from the method above
M_gradSpline = E*I.*spline_gradient_2';


end