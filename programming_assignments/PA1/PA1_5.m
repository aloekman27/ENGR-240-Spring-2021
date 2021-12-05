SumTerms = leibnizSummation(6);

function sumTerms = leibnizSummation(N)
%Enter the commands for your function here. Be sure to assign a value
%to the output variable sumTerms defined in the function command on line 1.

%n is the number of terms we want from the first term
n = 1:N;

%use build in sum to add the first term through the Nth term
sumTerms = sum(((-1).^(n-1))./((2*n)-1))

end