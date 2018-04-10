% Math 226B - HW #3
% Problem 5d
% Write a Matlab function that uses J, I, and VL from (b) to efficiently 
% compute the solution of upper-triangular linear systems with coefficient 
% matrix L^T .


function [x] = Usolve(J,I,VL,c)

n = length(c);
x = zeros(n,1);

x(n) = c(n)./VL(I(n));

for k = n-1:-1:1
    indexU = I(k):I(k+1) - 1;
    indexU = indexU(2:end);
    rowIndU = J(indexU);
    x(k) = c(k) - sum(VL(indexU).*x(rowIndU));
    x(k) = x(k)/VL(I(k));
end
end
    