% Math 226B - Homework #4
% Problem 3
% Use the built-in GMRES routine to write a Matlab program that lets you 
% run (full) GMRES or restarted GMRES (with restart parameter k0) for any 
% system (1), using the zero vector x0 = 0 as initial guess. As output, 
% your routine should produce the complete history of all the relative 
% residual norms produced during each run, the final approximate solution 
% xk of (1), and the total number of matrix-vector products q = Av computed 
% during each run.

function x = GMRESfunct(tol)

format long e
load('HW4_Problem3.mat')

n = length(b);
x0 = zeros(n,1);
maxit = n;
k0Vec = [2,5,10,20,50,100];

% full GMRES, no restart
[x,flag,relres,iter,resvec] = gmres(A,b,[],tol,maxit,[],[],x0);
iter
rel_resid = resvec./(revec(1));
tot_its = length(rel_resid)
k = (1:length(rel_resid));

hold on
plot(k,log(rel_resid),'LineWidth',1.1)
xlabel('Iteration Number')
ylabel('Log of Relative Residual')
    
% restarted GRMES   
for i=1:length(k0Vec)
    k0=k0Vec(i);
    [x,flag,relres,iter,resvec] = gmres(A,b,k0,tol,maxit,[],[],x0);
    iter
    rel_resid = resvec./(resvec(1));
    tot_its = length(rel_resid)
    k = (1:length(rel_resid));
    
    plot(k,log(rel_resid),'LineWidth',1.1)
    
end
legend('No Restart','k0 = 2','k0 = 5','k0 = 10','k0 = 20','k0 = 50','k0 = 100')
end
    
