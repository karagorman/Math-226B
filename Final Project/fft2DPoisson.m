
% Math 226B - Homework #3
% Problem 1b
% Write a Matlab program that implements your FFT-based solver from (a) 
% using Matlab?s ?fft?. The inputs should be the integer m determining 
% the mesh size h := 1/(m+1), an m×m matrix of the values fjk := f(xj,yk) 
% of the right-hand side function f at the grid points (xj,yk) := (jh,kh), 
% and 4 vectors of length m containing the values of the functions b0, b1, 
% c0, and c1 at the grid points on the 4 pieces of the boundary. The 
% output should be an m × m matrix of approximate solutions vjk.

function V = fft2DPoisson(m,f,c,ic,d,id)

format long e

h = 1/(2^m);

mx = (2^m)*ic - 1;
my = (2^m)*id - 1;

hx = 1/(mx + 1);
hy = 1/(my + 1);
alpha = 1;

% compute f'=z^T*f*z
f = fftMult(f); 
f = fftMult(f.').'; 

lambdax = 2*(1 - cos(pi*hx*(1:mx)));
lambday = 2*(1 - cos(pi*hy*(1:my)));

V_bar = zeros(mx,my);

for j = 1:mx
    V_bar(j,:) = f(j,:)./(lambdax(j) + alpha.*lambday);
end

% compute V=z*V'*z^T
V = fftMult(V_bar); 
V = fftMult(V.').'; 

% function to do matrix-vector multiplication using fft
    function w = fftMult(A)
        [m1,n] = size(A); hh = 1/(m1 + 1);
        A_tilde = [zeros(1,n); A; zeros(m1+1,n)];
        w_tilde = fft(A_tilde);
        w_hat = w_tilde(2:m1+1,:);
        w = -sqrt(2*hh)*imag(w_hat);
    end
end




