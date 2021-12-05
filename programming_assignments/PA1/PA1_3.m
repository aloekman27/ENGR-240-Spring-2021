[Up, UpDown] = vectorFun(5);

function [Up, UpDown] = vectorFun(A)
%Enter the commands for your function here. Be sure to assign values 
%to each of the output variables defined in the function command on line 1.

%create the first vector.
Up = [(0:2:2*A)] % 0 is the first element, increment is 2 (even numbers), 2*A is the last element

% this is the combination of the first vector(Up) with 2*A-1, -2 increments, and 1 is the last element
UpDown = [Up,(2*A-1:-2:1)]

end