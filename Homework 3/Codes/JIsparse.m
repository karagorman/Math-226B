% Math 226B - HW #3
% Problem 5a
% Write a Matlab function that efficiently generates the integer vectors 
% J and I for a given sparse matrix A.


function [J,I,VA] = JIsparse(A)

AL = tril(A);
[J,K,VA] = find(AL);

I(1) = 1;

for i = 2:size(AL,2)+1
   count = nnz(AL(:,i-1));
   I(i) = I(i-1) + count;
end
I = transpose(I);
end