% Math 226B - HW #3
% Problem 5e
% Use Matlab?s ?pcg? to compute the solution of symmetric positive 
% definite linear systems Ax = b using the CG method without 
% preconditioning and with the incomplete Cholesky preconditioner 
% generated by your function from (b). Employ your functions from 
% (c) and (d) for the solution of linear systems with L and LT .

function [J,I,VL,xNoP,xP,relresNoP,iterNoP,relresP,iterP] = SPDpcgSolver

% load('pcg_small.mat')
load('pcg_large.mat')

tol = 10^(-9);
n = length(b);
maxit = n;
x0 = ones(n,1);

% call function to generate J, I, and VA
[J,I,VA] = JIsparse(A);

% call function to perform incomplete Cholesky factorization
[VL] = incompleteCholesky(J,I,VA,A);

% PCG without preconditioning
[xNoP,flagNoP,relresNoP,iterNoP] = pcg(A,b,tol,maxit,[],[],x0);

% PCG with preconditioning
[xP,flagP,relresP,iterP] = pcg(A,b,tol,maxit,@(x) Lsolve(J,I,VL,x),...
    @(x) Usolve(J,I,VL,x),x0);

end

