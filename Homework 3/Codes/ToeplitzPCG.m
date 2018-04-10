% Math 226B - HW #3
% Problem 4a
% Write a Matlab routine based on your algorithm from Problem 5(d)of 
% Homework 1 and Matlab?s ?pcg? for solving linear systems (6) by means 
% of the CG method.


% woohoo! this works! (no idea if the answer is correct but it runs and
% converges quickly with a small residual)

function [x] = ToeplitzPCG(p,n)

format long e 

tol = 10^(-9);
maxit = n;

b = ones(n,1);

i=(1:n);
t = 1./((1 + sqrt(i-1)).^p);

[x,flag,relres,iter] = pcg(@TmultFunct,b,tol,maxit);
x1=x(1)
x2=x(100000)
x3=x(500000)
x4=x(700000)
x5=x(1000000)
relres
iter

    % function to do matrix-vector multiplication using FFT
    function y = TmultFunct(z)
        y = fftToeplitz(t,z);
        
    end
end