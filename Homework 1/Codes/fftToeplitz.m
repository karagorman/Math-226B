% Math 226B - Homework #1
% Problem 5d
% Use (c) to devise an algorithm that computes matrix-vector products 
% y = T x with Toeplitz matrices T via your Matlab function form (b).
% Write a Matlab function that implements this approach.

function [y] = fftToeplitz(t,x)

n = length(x);

c = zeros(1,2*n-1);

c(1:n) = t(n:2*n-1);
c(n+1:2*n-1) = t(1:n-1);

x_tilde = zeros(2*n-1,1);

for i = 1:n
    x_tilde(i) = x(i);
end

c = c';

y_tilde = fftCirculant(c,x_tilde);

y = y_tilde(1:n);

end


 
