A = [3 17 4 5 10 15];
B = [33 12 6 31 37 27 49 22 13 28]';
[ABrow, BAcol, FirstHalfA_LastHalfB] = vectorFun(A, B)

function [ABrow, BAcol, FirstHalfA_LastHalfB] = vectorFun(A, B)
%Enter the commands for your function below Be sure to assign the values to each of the output variables.
%Defined in the function command on line 1.

%create the vector consisting from vectors A and B
ABrow = [A,B']
%use ABrow and flip its order and change to column vectors from row vectors
BAcol = [(fliplr(ABrow))']
%Use the first three elements from A and the last 5 elements from B to generate a new row vector
FirstHalfA_LastHalfB = [A(1:3),B(6:10)']

end