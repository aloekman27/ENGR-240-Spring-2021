[terms, seriesSums] = maclaurinCosine(pi/3, 5)

function [terms, seriesSums, approxRelError, truRelError] = maclaurinCosine(x, N)
%Enter the commands for your function here.  
 

% set up a column vector from the first term, 0, to the Nth term. The first
% term is 0 because it effects the exponent and factorial
N = (0:N-1)';
terms = ((-1).^N .* x.^(2.*N)) ./ factorial(2.*N); %use the formula to calculate each term of the serial to the Nth term

%sum up the values from the first term up to the fifth term and put their
%sums in order in a column vector
seriesSums = cumsum(terms);

%
%calculate the approximate relative error of the series
presApprox = (seriesSums(2:end)); %calculate the present approximation(the 2nd term of the series until the end of the series)
prevApprox = (seriesSums(1:end-1)); %calculate the previous approximation(the first term up to the 2nd to last term) 

%using the formula to find the approximate relative error of the
%calculation
approxRelError = abs((presApprox - prevApprox)./(presApprox));

%
%calculate the true relative error of the series
trueValue = cos(x); %find the theoretical value of the cos(x)
approx = seriesSums(1:end); %the terms of the maclaurin series

%use the equation below to compare the true relative error using the true
%value and the approximation. as we have more terms the true relative error
%value decreases and gets close to 0
truRelError = abs((trueValue-approx)/trueValue);

end