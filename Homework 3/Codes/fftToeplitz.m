% Math 226B - Homework #1
% Problem 5d
% Use (c) to devise an algorithm that computes matrix-vector products 
% y = T x with Toeplitz matrices T via your Matlab function form (b).
% Write a Matlab function that implements this approach.

function [y] = fftToeplitz(t,x)

m = length(x);

c = zeros(1,2*m-1);

c(1:m) = t(1:m);
c(m+1:2*m-1) = t(m:-1:2);

x_tilde = zeros(2*m-1,1);

x_tilde(1:m) = x;

%c = c';

y_tilde = fftCirculant(c,x_tilde);


y = y_tilde(1:m);

end


 
