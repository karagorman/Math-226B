% Math 226B - HW #3
% Problem 4b
% Based on Matlab's pcg, write a Matlab routine for solving linear 
% systems C_T^-1x = C_T^-1b by means of the preconditioned CG method.

function [x,relres,iter] = ToepCircPCG(p,n)

format long e

tol = 10^(-9);
maxit = n;
b = ones(n,1);

i=(0:n-1);
t = 1./((1 + sqrt(i)).^p);

c = ((n:-1:1).*t + (0:(n-1)).*t([1,end:-1:2]))./n;


[x,flag,relres,iter] = pcg(@(x) TmultFunct(t,x),b,tol,maxit,@(x) cTimult(c,x))

% x1=x(1)
% x2=x(100000)
% x3=x(500000)
% x4=x(700000)
% x5=x(1000000)
% relres=relres
% iter=iter

    % function to multiply matrix C_T^-1 and vector v
    function h = cTimult(c,v)
        lambdaVec = 1./conj(fft(c'));
        h = conj(fft(conj(lambdaVec.*fft(v))))/n;
    end

    % function to multiply topelitz matrix T and z
    function y = TmultFunct(t,z)
        y = fftToeplitz(t,z);
    end        
end