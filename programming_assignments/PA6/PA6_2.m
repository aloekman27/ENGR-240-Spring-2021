%Example function call for first four tests
xDat = linspace(-pi/2,pi/2,7);
[yjPoly, yjLinear, yjClamped, Et_matrix]=PA6_interp(xDat, 0.02)

function [yjPoly, yjLinear, yjClamped, Et_matrix]=PA6_interp(xSamples, increment)
% [yjPoly, yjLinear, yjClamped, Et_matrix]=PA6_interp(xSamples, increment)
%
% Interpolation Comparison where 3 types of interpolation approaches is used on a 
% known mathematical function and compared with absolute errors  corresponding
% to each point in the graph
%
% Inputs: 
%   xSamples  : A vector of x values to generate the sampled data
%   increment : A scalar specifying the increment that will be used in interpolation
%
% Outputs:
%   yjPoly    : The interpolated data generated with polynomial interpolation
%   yjLinear  : The interpolated data generated with a linear spline
%   jyClamped : The interpolated data generated with the clamped cubic spline
%   Et_matrix : A 3-column matrix of the absolute error for each approach

%Define the variables that will be the x, fx, and interpolated xx values
xSamples'; %x values
fx = exp(-xSamples.^2).*cos(2.*xSamples); %f(x) values
xx = min(xSamples):increment:max(xSamples);
n=length(xSamples)-1; %define the N-1 polynomial order where n is the points in the sampled data

%find the f(x) from the xx values which will create a smooth line in the graph
fx_error= (exp(-xx.^2).*cos(2.*xx))';

%Polynomial Interpolation
yjPoly_Coeffs = polyfit(xSamples,fx,n);
yjPoly = polyval(yjPoly_Coeffs,xx)';

%Linear Spline Interpolation
yjLinear = interp1(xSamples,fx,xx,'linear')';

%Cubic Spline with Clamped Endpoints
%derivatives to the function f(x) evaluated at  the left end point
d_fx_min = -2*exp(-min(xSamples)^2)*min(xSamples)*cos(2*min(xSamples)) - 2*exp(-(min(xSamples))^2)*sin(2*min(xSamples));
%derivatives to the function f(x) evaluated at  the right end point
d_fx_max = -2*exp(-max(xSamples)^2)*max(xSamples)*cos(2*max(xSamples)) - 2*exp(-(max(xSamples))^2)*sin(2*max(xSamples));
fx_clamp = [d_fx_min, fx, d_fx_max]; %row vector of clamped conditions
yjClamped = spline(xSamples,fx_clamp,xx)';

%Calculate the absolute error for each point of the approaches used
Et_yjPoly = (fx_error-yjPoly)
Et_yjLinear = (fx_error-yjLinear)
Et_yjClamped = (fx_error-yjClamped)

%Allocate the absolute error values into a matrix
Et_matrix = [Et_yjPoly Et_yjLinear Et_yjClamped]

%Plot the graph of the function and interpolation
figure(1) 
plot(xSamples,fx,'bd')
hold on
plot(xx,yjPoly,'k:') %6th order polynomial graph
plot(xx,yjLinear,'r:') %piecewise linear spline graph
plot(xx,yjClamped,'g:')%spline with clamped endpoints graph
plot(xx,fx_error,'m--')%smooth curve of data
title('Interpolation of a Known Mathematical Function')
xlabel('x')
ylabel('f(x)')
legend('data', '6^{th} Order Polynomial','Piecewise Linear','Spline (clamped)','e^{-xx{^2}}*cos(2*xx))')
grid
hold off

end