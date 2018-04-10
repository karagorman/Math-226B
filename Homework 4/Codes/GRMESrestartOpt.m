% Math 226B _ Homework #4
% Problem 5
% Use Matlab?s function ?gmres? to write Matlab programs for each of the 
% following algorithms for solving linear systems (1):
% (i) GMRES (without preconditioning);
% (ii) Restarted GMRES (without preconditioning);
% (iii) GMRES with diagonal preconditioning applied from the right, 
% i.e., M1 =I and M2 =D0, where D0 denotes the diagonal part of A, see (3);
% (iv) Restarted GMRES with diagonal preconditioning applied from the right;
% (v) GMRES with the SSOR-type preconditioner derived in Problem 4;
% (vi) Restarted GMRES with the SSOR-type preconditioner derived in 
% Problem 4.

function x = GRMESrestartOpt(fileNum)

if (fileNum == 1)
    load('HW4_Problem5b_1.mat')
elseif (fileNum == 2)
    load('HW4_Problem5b_2.mat')
end

n = length(b);
x0 = ones(n,1);
maxit = n;
tol = 1e-8;
k0Vec = [5,10,20];


% Case 1: GMRES without preconditioning
[x,flag,relres,iter,resvec] = gmres(A,b,[],tol,maxit,[],[],x0);
iter
rel_resid = resvec./resvec(1);
tot_its = length(rel_resid)
k = (1:length(rel_resid));

hold on
subplot(4,2,1)
plot(k,log(rel_resid),'LineWidth',1)
xlabel('Iteration Number')
ylabel('Log of Relative Residual')
title('Full GMRES, No Preconditioning')

    
% Case 2: Restarted GMRES without preconditioning
k0Vec = [5,10,20];
for i=1:length(k0Vec)
    k0=k0Vec(i);
    [x,flag,relres,iter,resvec] = gmres(A,b,k0,tol,maxit,[],[],x0);
    iter
    rel_resid = resvec./resvec(1);
    tot_its = length(rel_resid)
    k = (1:length(rel_resid));

    subplot(4,2,2)
    hold on
    plot(k,log(rel_resid),'LineWidth',1)
    xlabel('Iteration Number')
    ylabel('Log of Relative Residual')
    title('Restarted GMRES, No Preconditioning')
end
hold off
legend('k0 = 5','k0 = 10','k0 = 20')

% Case 3: GMRES with diagonal preconditioning applied from the right
D0 = diag(diag(A));
I = speye(n);
[x,flag,relres,iter,resvec] = gmres(A,b,[],tol,maxit,I,D0,x0);
iter
rel_resid = resvec./resvec(1);
tot_its = length(rel_resid)
k = (1:length(rel_resid));

subplot(4,2,3)
plot(k,log(rel_resid),'LineWidth',1)
xlabel('Iteration Number')
ylabel('Log of Relative Residual')
title('Full GMRES, With Diagonal Preconditioning')
    
% Case 4: Restarted GMRES with diagonal preconditioning applied from the
% right
for i=1:length(k0Vec)
    k0=k0Vec(i);
    [x,flag,relres,iter,resvec] = gmres(A,b,k0,tol,maxit,I,D0,x0);
    iter
    rel_resid = resvec./resvec(1);
    tot_its = length(rel_resid)
    k = (1:length(rel_resid));

    subplot(4,2,4)
    hold on
    plot(k,log(rel_resid),'LineWidth',1)
    xlabel('Iteration Number')
    ylabel('Log of Relative Residual')
    title('Restarted GMRES, With Diagonal Preconditioning')
end
hold off
legend('k0 = 5','k0 = 10','k0 = 20')

% Case 5: Full GMRES with SSOR-type preconditioning (from problem 4)
D0 = diag(diag(A));
F = -tril(A,-1);
G = -triu(A,1);
I = speye(n);

%D=D0;
D = D0;
maxit = n;
D1 = D0 - 2*D;
L = D-F;
U = D-G;
M1 = L*D^(-1);
M2 = U;
bp = M1\b;
x0p = M2*x0;
% if we want the actual solution x_k, need to multiply it on the right by
% M2^-1 at the end

[x,flag,relres,iter,resvec] = gmres(@(v) ApMultFunct(L,U,D,D1,v),bp,[],tol,maxit,[],[],x0p);
iter
rel_resid = resvec./resvec(1);
tot_its = length(rel_resid)
k = (1:length(rel_resid));

subplot(4,2,5)
plot(k,log(rel_resid),'LineWidth',1)
xlabel('Iteration Number')
ylabel('Log of Relative Residual')
title('Full GMRES, With SSOR-type Preconditioning, with D=D0')   

% D=10I
D = 10*I;
D1 = D0 - 2*D;
L = D-F;
U = D-G;
M1 = L*D^(-1);
M2 = U;
bp = M1\b;
x0p = M2*x0;
    
[x,flag,relres,iter,resvec] = gmres(@(v) ApMultFunct(L,U,D,D1,v),bp,[],tol,maxit,[],[],x0p);
iter
rel_resid = resvec./resvec(1);
tot_its = length(rel_resid)
k = (1:length(rel_resid));

subplot(4,2,6)
plot(k,log(rel_resid),'LineWidth',1)
xlabel('Iteration Number')
ylabel('Log of Relative Residual')
title('Full GMRES, With SSOR-type Preconditioning, with D=10I')

% Case 6: Restarted GMRES with SSOR-type preconditioning
%D=D0
D = D0;
maxit = n;
D1 = D0 - 2*D;
L = D-F;
U = D-G;
M1 = L*D^(-1);
M2 = U;
bp = M1\b;
x0p = M2*x0;

for i=1:length(k0Vec)
    k0=k0Vec(i); 

    [x,flag,relres,iter,resvec] = gmres(@(v) ApMultFunct(L,U,D,D1,v),bp,k0,tol,maxit,[],[],x0p);
    iter
    rel_resid = resvec./resvec(1);
    tot_its = length(rel_resid)
    k = (1:length(rel_resid));

    subplot(4,2,7)
    hold on
    plot(k,log(rel_resid),'LineWidth',1)
    xlabel('Iteration Number')
    ylabel('Log of Relative Residual')
    title('Restarted GMRES, With SSOR-type Preconditioning, with D = D0')
end
hold off
legend('k0 = 5','k0 = 10','k0 = 20')
 
%D=10I
D = 10*I;
D1 = D0 - 2*D;
L = D-F;
U = D-G;

M1 = L*D^(-1);
M2 = U;
bp = M1\b;
x0p = M2*x0;

for i=1:length(k0Vec)
    k0=k0Vec(i); 

    [x,flag,relres,iter,resvec] = gmres(@(v) ApMultFunct(L,U,D,D1,v),bp,k0,tol,maxit,[],[],x0p);
    iter
    rel_resid = resvec./resvec(1);
    tot_its = length(rel_resid)
    k = (1:length(rel_resid));

    subplot(4,2,8)
    hold on
    plot(k,log(rel_resid),'LineWidth',1)
    xlabel('Iteration Number')
    ylabel('Log of Relative Residual')
    title('Restarted GMRES, With SSOR-type Preconditioning, with D = 10I')
end
hold off
legend('k0 = 5','k0 = 10','k0 = 20')

   
end
    