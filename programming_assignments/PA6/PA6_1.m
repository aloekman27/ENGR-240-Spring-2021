%Example function call for first two tests 
[erfLagrange2, EtLagrange2, erfSpline, EtSpline]  = PA6_erf(0.92)

function [erfLagrange2, EtLagrange2, erfSpline, EtSpline] = PA6_erf(x)
% [erfLagrange2, EtLagrange2, erfSpline, EtSpline] = PA6_erf(x)
%
% Use 2nd Order Lagrange Polynomial and Spline "not-a-knot" Interpolation
% to evaluate an error function from a tabulated data
%
% Input :
%   x : scalar input that will be evaluated using interpolation
%   
% Output : 
%   erfLagrange2 : An interpolated value erx(x) using a 2nd order Lagrange
%                  interpolating polynomial
%   EtLagrange2  : The true error associated with erfLangrange2
%   erfSpline    : An interpolated value erx(x) using a Spline "not-a-knot"
%                  end conditions
%   EtSpline     : The true error associated with erfSpline

%Define values for initial x and erfx corresponding to each x
xi = [0 0.4 0.8 1.2 1.6 2];
fx = [0 0.42839 0.7421 0.91031 0.97635 0.99532];

%These if else statements are used for the 2nd Order Lagrange Polynomial

if x >= xi(1) && x <= xi(2) %if x is between 0 and 0.4
    %define coefficients
    x1 = xi(1);
    x2 = xi(2);
    x3 = xi(3);
    %Calculate the weighting coefficients
    L1 = ((x-x2)*(x-x3))/((x1-x2)*(x1-x3));
    L2 = ((x-x1)*(x-x3))/((x2-x1)*(x2-x3));
    L3 = ((x-x2)*(x-x1))/((x3-x2)*(x3-x1));
    
    %Calcualte the Lagrange interpolation
    erfLagrange2 = L1*fx(1)+L2*fx(2)+L3*fx(3);
    %Calculate the true error associated with the Lagrange Interpolation
    EtLagrange2 = erf(x) - erfLagrange2; 
    
        %if x is between 0.4 and 0.8, third xi is the closest to data x
elseif x >= xi(2) && x <= xi(3) & (xi(1)-x < xi(4)-x) 
    x1 = xi(2);
    x2 = xi(3);
    x3 = xi(1);
    
    L1 = ((x-x2)*(x-x3))/((x1-x2)*(x1-x3));
    L2 = ((x-x1)*(x-x3))/((x2-x1)*(x2-x3));
    L3 = ((x-x2)*(x-x1))/((x3-x2)*(x3-x1));
    
    
    erfLagrange2 = L1*fx(2)+L2*fx(3)+L3*fx(1);

    EtLagrange2 = erf(x) - erfLagrange2;
    
    %if x is between 0.4 and 0.8, third xi is the closest to data x
elseif x >= xi(2) && x <= xi(3) & (xi(1)-x > xi(4)-x) 
    x1 = xi(2);
    x2 = xi(3);
    x3 = xi(4);
    
    L1 = ((x-x2)*(x-x3))/((x1-x2)*(x1-x3));
    L2 = ((x-x1)*(x-x3))/((x2-x1)*(x2-x3));
    L3 = ((x-x2)*(x-x1))/((x3-x2)*(x3-x1));
    
    
    erfLagrange2 = L1*fx(2)+L2*fx(3)+L3*fx(4);

    EtLagrange2 = erf(x) - erfLagrange2;

    %if x is between 0.8 and 1.2, third xi is the closest to data x
elseif x >= xi(3) && x <= xi(4) & (xi(2)-x < xi(5)-x)
    x1 = xi(3);
    x2 = xi(4);
    x3 = xi(2);
    
    L1 = ((x-x2)*(x-x3))/((x1-x2)*(x1-x3));
    L2 = ((x-x1)*(x-x3))/((x2-x1)*(x2-x3));
    L3 = ((x-x2)*(x-x1))/((x3-x2)*(x3-x1));
    
    
    erfLagrange2 = L1*fx(3)+L2*fx(4)+L3*fx(2);

    EtLagrange2 = erf(x) - erfLagrange2;
    
    %if x is between 0.8 and 1.2, third xi is the closest to data x
elseif x >= xi(3) && x <= xi(4) & (xi(2)-x > xi(5)-x)
    x1 = xi(3);
    x2 = xi(4);
    x3 = xi(5);
    
    L1 = ((x-x2)*(x-x3))/((x1-x2)*(x1-x3));
    L2 = ((x-x1)*(x-x3))/((x2-x1)*(x2-x3));
    L3 = ((x-x2)*(x-x1))/((x3-x2)*(x3-x1));
    
    
    erfLagrange2 = L1*fx(3)+L2*fx(4)+L3*fx(5);

    EtLagrange2 = erf(x) - erfLagrange2;  

    %if x is between 1.2 and 1.6, third xi is the closest to data x
elseif x >= xi(4) && x <= xi(5) & (xi(3)-x < xi(6)-x)
    x1 = xi(4);
    x2 = xi(5);
    x3 = xi(3);
    
    L1 = ((x-x2)*(x-x3))/((x1-x2)*(x1-x3));
    L2 = ((x-x1)*(x-x3))/((x2-x1)*(x2-x3));
    L3 = ((x-x2)*(x-x1))/((x3-x2)*(x3-x1));
    
    
    erfLagrange2 = L1*fx(4)+L2*fx(5)+L3*fx(3);

    EtLagrange2 = erf(x) - erfLagrange2;  
 
    %if x is between 1.2 and 1.6, third xi is the closest to data x
elseif x >= xi(4) && x <= xi(5) & (xi(3)-x > xi(6)-x)
    x1 = xi(4);
    x2 = xi(5);
    x3 = xi(6);
    
    L1 = ((x-x2)*(x-x3))/((x1-x2)*(x1-x3));
    L2 = ((x-x1)*(x-x3))/((x2-x1)*(x2-x3));
    L3 = ((x-x2)*(x-x1))/((x3-x2)*(x3-x1));
    
    
    erfLagrange2 = L1*fx(4)+L2*fx(5)+L3*fx(6);

    EtLagrange2 = erf(x) - erfLagrange2; 
    
    %if x is between 1.6 and 2, third xi is the closest to data x
else
    
    x1 = xi(5);
    x2 = xi(6);
    x3 = xi(4);
    
    L1 = ((x-x2)*(x-x3))/((x1-x2)*(x1-x3));
    L2 = ((x-x1)*(x-x3))/((x2-x1)*(x2-x3));
    L3 = ((x-x2)*(x-x1))/((x3-x2)*(x3-x1));
    
    
    erfLagrange2 = L1*fx(5)+L2*fx(6)+L3*fx(4);

    EtLagrange2 = erf(x) - erfLagrange2; 
        
end

%The code below calculates interpolation using splines with
%"not-a-knot" end condition

%calculate spline interpolation at point x
erfSpline = spline(xi,fx,x);
EtSpline = erf(x) - erfSpline; %error associated with the spline interpolation



end