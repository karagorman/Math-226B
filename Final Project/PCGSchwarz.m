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


function [v_pcgs,V_pcgs,u,NumSolvesR1,NumSolvesR2] = PCGSchwarz(tol,nmax,...
    m,c,d,alpha,beta,gamma)

format long e

[A,b,A1,A2,I1,I2,Num_Mat,tot_R1_pts,tot_R2_pts,tot_pts,R1_xind,R1_yind,...
    R2_xind,R2_yind] = AbData(m,c,d,alpha,beta,gamma);

v0 = zeros(length(b),1);

h = 1/(2^m);
n = length(b);

v0p = v0;
% compute b'=M^-1*b
s1 = I1*b;
s2 = I2*b;
kmin1 = 1;
kmax1 = tot_R1_pts;
w1 = RiSolve(R1_xind,R1_yind,kmin1,kmax1,s1,0,1,0,3,Num_Mat,m);
NumSolvesR1 = 1;
kmin2 = tot_pts - tot_R2_pts + 1;
kmax2 = tot_pts;
w2 = RiSolve(R2_xind,R2_yind,kmin2,kmax2,s2,0,4,0,1,Num_Mat,m);
NumSolvesR2 = 1;
bp = (I1.')*w1 + (I2.')*w2;


v_pcgs = pcg(@(v) ApMult(A,v),bp,tol,nmax,[],[],v0p);

V_pcgs = zeros(size(Num_Mat,1),size(Num_Mat,2));

for i = 1:n
    [data_x,data_y] = find(Num_Mat == i);
    V_pcgs(data_x,data_y) = v_pcgs(i); 
end

x=0:h:c+4;
y=0:h:3;
[X,Y] = meshgrid(x,y);
X= X'; Y = Y';
u = zeros(size(Num_Mat));
for i=1:n
   [data_x,data_y] = find(Num_Mat == i);
   u(data_x,data_y) = Y(data_x,data_y).^(alpha).*...
       sin(beta*pi*X(data_x,data_y)).*cos(gamma*pi*Y(data_x,data_y));
end

% function to compute q = A'p
    function q = ApMult(A,v)
        v_bar = A*v;
        s1 = I1*v_bar;
        s2 = I2*v_bar;
        kmin1 = 1;
        kmax1 = tot_R1_pts;
        w1 = RiSolve(R1_xind,R1_yind,kmin1,kmax1,s1,0,1,0,3,Num_Mat,m);
        kmin2 = tot_pts - tot_R2_pts + 1;
        kmax2 = tot_pts;
        w2 = RiSolve(R2_xind,R2_yind,kmin2,kmax2,s2,0,4,0,1,Num_Mat,m);
        q = (I1.')*w1 + (I2.')*w2;
        NumSolvesR1 = NumSolvesR1 + 1;
        NumSolvesR2 = NumSolvesR2 + 1;
    end
end