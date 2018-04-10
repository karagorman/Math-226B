% Math 226B - Final Project
% Part 3b
% For each of the following two domain decomposition methods, write a 
% Matlab program that solves the discretized problem Av = b:
% CG method with additive Schwarz preconditioning.
% The program for each of these methods should allow as inputs an arbitrary 
% initial guess v(0) for the solution of Av = b, a ?small? convergence 
% tolerance tol for the domain decomposition method, and a large integer 
% nmax to be used as a safeguard against an excessive number of iterations 
% of the domain decomposition method.
% For method (2), apply the preconditioner M from the left, i.e., 
% solve M?1Av = M?1b, and employ Matlab?s ?pcg? function.


function [v_pcgs,V_pcsg] = PCGSchwarz(v0,tol,nmax,m,A,A1,A2,b,I1,I2)

[A,b,A1,A2,I1,I2,Num_Mat,tot_R1_pts,tot_R2_pts,tot_pts] = AbData(m,c,d,alpha,beta,gamma);

h = 1/(2^m);
n = length(b);


% use fft solver for part of these?
B1 = (I1.')*(A1\I1);
B2 = (I2.')*(A2\I2);



Ap = (B1 + B2)*A;
bp = (B1 + B2)*b;
v0p = v0;


[v_pcgs,flag,relres,iter,resvec] = pcg(Ap,bp,[],[],v0p);

V_pcsg = zeros(size(Num_Mat,1),size(Num_Mat,2));

for i = 1:n
    [data_x,data_y] = find(Num_Mat == i);
    V_pcsg(data_x,data_y) = v_pcsg(i); 
end






% find relative error
% [X,Y]=meshgrid(0:h:c+4,0:h:3);
% X = X'; Y = Y';
% u = Y.^(alpha).*sin(beta*pi*X).*cos(gamma*pi*Y);
% 
% [max_V, max_ind] = max(V_pcgs(:));
% [row_ind,col_ind] = ind2sub(size(V_pcgs),max_ind);
% rel_error_pcgs = abs(max_V - u(row_ind,col_ind))/abs(u(row_ind,col_ind));


end