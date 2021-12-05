%This script m-file calculates the summation above 

%Set up values that will make up the equation from the information above
k = 5; m = 500000; n = 0.00001;



%calculate the value of x1, this time we will use double precision
x1 = (k-(m*n));

%set up the variable ZeroDouble to call the value of x using double precision
ZeroDouble=x1-0


%%set variables and set them as the single precision values of each values of k, m, and n
a = single(k);
b = single(m);
c = single(n);

%use repmat to do the summation of the values in the equation
d = repmat(c,b,1);
%set the value of the summation to single precision
e = single(sum(d));

%calculate the equation above
sum = a-e;
%set up the variable ZeroDouble to call the value of x using single precision
ZeroSingle = single(sum)