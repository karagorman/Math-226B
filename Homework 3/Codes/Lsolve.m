% Math 226B - HW #3
% Problem 5c
% Write a Matlab function that uses J, I, and VL from (b) to efficiently 
% compute the solution of lower-triangular linear systems with coefficient 
% matrix L.

function [c] = Lsolve(J,I,VL,b)

n = length(b);
c = zeros(n,1);

for k = 1:n-1
    c(k) = b(k);
    indexL = I(k):I(k+1) - 1;
    rowIndL = J(indexL);
    b(rowIndL) = b(rowIndL) - (VL(indexL)*c(k))./VL(indexL(1));
end
c(n) = b(n)/VL(end);
end
