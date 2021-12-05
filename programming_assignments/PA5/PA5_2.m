%load file with the data.  Note this file is available on Canvas
load PA5thermis.dat 
[abc, rSquared, standardError] = PA5_thermistor(PA5thermis)

function [abc, rSquared, standardError] = PA5_thermistor(data)
%[abc, rSquared, standardError] = PA5_thermistor(data)
%
% A function that uses General Linear Least Squares Approach to determine
% the best values of the constants a, b, and c to fit the model
%
%Inputs: 
%   data            : A two-column matrix with resistance values in the first column
%                 and temperature values in the second column
%
%Outputs:
%   abc             : A 1x3 row vector of the fit coefficients a, b, and c
%   rSquared        : The coefficient of determination for the curve fit
%   standard Error  : The standard error for the fit


%Define variables for our inputs, resistance and temperature
R = data(:,1)'; %Ohms
T = data(:,2)'; %Celcius
y=1./T; %linearized y



%Basis Functions
Z0 = @(x) log(R).^0;
Z1 = @(x) log(R);
Z2 = @(x) log(R).^3;
Z = [Z0(R)', Z1(R)', Z2(R)'];

%Solve for coefficients a, b, and c
A = Z\y';
abc = A';


%Set up a function to solve for T from the coefficients we have found
T_model = @(R) 1./(abc(1).*log(R).^0 + abc(2).*log(R).^1 + abc(3).*log(R).^3);

St = sum((T-mean(T)).^2);
Sr = sum((T-T_model(R)).^2);%sum of the squared residuals
rSquared = 1-Sr/St; %the sum of the squares of deviation from the mean


%Calculate the standard error of fit
n = length(R); %the data points provided from the data
d = length(abc); %the degrees of freedom (ie. the number of coefficients)
standardError = sqrt(Sr/(n-d)); 

end