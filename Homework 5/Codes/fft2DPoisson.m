
% Math 226B - Homework #3
% Problem 1b
% Write a Matlab program that implements your FFT-based solver from (a) 
% using Matlab?s ?fft?. The inputs should be the integer m determining 
% the mesh size h := 1/(m+1), an m×m matrix of the values fjk := f(xj,yk) 
% of the right-hand side function f at the grid points (xj,yk) := (jh,kh), 
% and 4 vectors of length m containing the values of the functions b0, b1, 
% c0, and c1 at the grid points on the 4 pieces of the boundary. The 
% output should be an m × m matrix of approximate solutions vjk.

function V = fft2DPoisson(m,f)

format long e

h = 1/(m+1);

% compute f'=z^T*f*z
f = fftMult(f);
f = fftMult(f.').';

[X,Y]=meshgrid(1:m,1:m);

lambda = 2*(1 - cos(pi*h.*X) + 1 - cos(pi*h.*Y));
V_bar = f./lambda;

% compute V=z*V'*z^T
V = fftMult(V_bar);
V = fftMult(V.').';

% function to do matrix-vector multiplication using fft
    function w = fftMult(A)
        
        n = size(A,2);
        A_tilde = [zeros(1,n); A; zeros(n+1,n)];
        w_tilde = fft(A_tilde);
        w_hat = w_tilde(2:n+1,:);
        w = -sqrt(2*h)*imag(w_hat);
        
    end
end




