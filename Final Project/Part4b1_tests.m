% To Run Part 4b1
clear all

%%
% case (i):
m = 6; c=2/5; d=1/5; alpha=0; beta=1; gamma=1/2;
h = 1/(2^m);
tol = 1e-9;
nmax = 100;

[u_exact,v_das,V_das1,V_das2,numR1solves,numR2solves,u1,u2] = ...
    DiscAltSchwarz(tol,nmax,m,c,d,alpha,beta,gamma);
numR1solves
numR2solves
erroridas = norm(v_das - u_exact,inf)/norm(u_exact,inf)

% plot disc alt schwarz
x=h*(1:size(V_das1,1));
y=h*(1:size(V_das1,2));
[X,Y] = meshgrid(x,y);
X= X'; Y = Y';

subplot(3,4,1)
hold on
surf(X,Y,V_das1,'EdgeColor','none')
surf(X,Y,V_das2,'EdgeColor','none')
axis tight
colorbar
title('Case (i): Disc. Alternating Schwarz')
hold off

[v_gmres,V_gmres1,V_gmres2,its,numR1solves,numR2solves] = ...
    GMRESAltSchwarz(tol,nmax,m,c,d,alpha,beta,gamma);
numR1solves
numR2solves
errorigmres = norm(v_gmres - u_exact,inf)/norm(u_exact,inf)

subplot(3,4,2)
hold on
surf(X,Y,V_gmres1,'EdgeColor','none')
surf(X,Y,V_gmres2,'EdgeColor','none')
axis tight
colorbar
title('Case (i): GMRES w/ Alt. Schwarz')
hold off





%%

%%
% case (ii)
m=6; c=2/5; d=1/5; alpha=2; beta=7/2; gamma=2;
tol = 1e-9;
nmax = 100;

[u_exact,v_das,V_das1,V_das2,numR1solves,numR2solves,u1,u2] = DiscAltSchwarz(tol,nmax,m,c,d,alpha,beta,gamma);
numR1solves
numR2solves
erroriidas = norm(v_das - u_exact,inf)/norm(u_exact,inf)

x=h*(1:size(V_das1,1));
y=h*(1:size(V_das1,2));
[X,Y] = meshgrid(x,y);
X= X'; Y = Y';

subplot(3,4,3)
hold on
surf(X,Y,V_das1,'EdgeColor','none')
surf(X,Y,V_das2,'EdgeColor','none')
axis tight
colorbar
title('Case (ii): Disc. Alternating Schwarz')
hold off

[v_gmres,V_gmres1,V_gmres2,its,numR1solves,numR2solves] = GMRESAltSchwarz(tol,nmax,m,c,d,alpha,beta,gamma);
numR1solves
numR2solves
erroriigmres = norm(v_gmres - u_exact,inf)/norm(u_exact,inf)

subplot(3,4,4)
hold on
surf(X,Y,V_gmres1,'EdgeColor','none')
surf(X,Y,V_gmres2,'EdgeColor','none')
axis tight
colorbar
title('Case (ii): GMRES w/ Alt. Schwarz')
hold off

%%

%%
% case (iii)
m=6; c=1/5; d=6/5; alpha=0; beta=1; gamma=1/2;

tol = 1e-9;
nmax = 100;

[u_exact,v_das,V_das1,V_das2,numR1solves,numR2solves,u1,u2] = DiscAltSchwarz(tol,nmax,m,c,d,alpha,beta,gamma);
numR1solves
numR2solves
erroriiidas = norm(v_das - u_exact,inf)/norm(u_exact,inf)

x=h*(1:size(V_das1,1));
y=h*(1:size(V_das1,2));
[X,Y] = meshgrid(x,y);
X= X'; Y = Y';

subplot(3,4,5)
hold on
surf(X,Y,V_das1,'EdgeColor','none')
surf(X,Y,V_das2,'EdgeColor','none')
axis tight
colorbar
title('Case (iii): Disc. Alternating Schwarz')
hold off

[v_gmres,V_gmres1,V_gmres2,its,numR1solves,numR2solves] = GMRESAltSchwarz(tol,nmax,m,c,d,alpha,beta,gamma);
numR1solves
numR2solves
erroriiigmres = norm(v_gmres - u_exact,inf)/norm(u_exact,inf)

subplot(3,4,6)
hold on
surf(X,Y,V_gmres1,'EdgeColor','none')
surf(X,Y,V_gmres2,'EdgeColor','none')
axis tight
colorbar
title('Case (iii): GMRES w/ Alt. Schwarz')
hold off



%%

%%
% case (iv)
m=6; c=1/5; d=6/5; alpha=2; beta=7/2; gamma=2;
tol = 1e-9;
nmax = 100;

[u_exact,v_das,V_das1,V_das2,numR1solves,numR2solves,u1,u2] = DiscAltSchwarz(tol,nmax,m,c,d,alpha,beta,gamma);
numR1solves
numR2solves
errorivdas = norm(v_das - u_exact,inf)/norm(u_exact,inf)

x=h*(1:size(V_das1,1));
y=h*(1:size(V_das1,2));
[X,Y] = meshgrid(x,y);
X= X'; Y = Y';

subplot(3,4,7)
hold on
surf(X,Y,V_das1,'EdgeColor','none')
surf(X,Y,V_das2,'EdgeColor','none')
axis tight
colorbar
title('Case (iv): Disc. Alternating Schwarz')
hold off

[v_gmres,V_gmres1,V_gmres2,its,numR1solves,numR2solves] = GMRESAltSchwarz(tol,nmax,m,c,d,alpha,beta,gamma);
numR1solves
numR2solves
errorivgmres = norm(v_gmres - u_exact,inf)/norm(u_exact,inf)

subplot(3,4,8)
hold on
surf(X,Y,V_gmres1,'EdgeColor','none')
surf(X,Y,V_gmres2,'EdgeColor','none')
axis tight
colorbar
title('Case (iv): GMRES w/ Alt. Schwarz')
hold off

%%

%%
% case (v)
m=6; c=2/5; d=9/5; alpha=0; beta=1; gamm=1/2;
tol = 1e-9;
nmax = 100;

[u_exact,v_das,V_das1,V_das2,numR1solves,numR2solves,u1,u2] = DiscAltSchwarz(tol,nmax,m,c,d,alpha,beta,gamma);
numR1solves
numR2solves
errorvdas = norm(v_das - u_exact,inf)/norm(u_exact,inf)

x=h*(1:size(V_das1,1));
y=h*(1:size(V_das1,2));
[X,Y] = meshgrid(x,y);
X= X'; Y = Y';

subplot(3,4,9)
hold on
surf(X,Y,V_das1,'EdgeColor','none')
surf(X,Y,V_das2,'EdgeColor','none')
axis tight
colorbar
title('Case (v): Disc. Alternating Schwarz')
hold off

[v_gmres,V_gmres1,V_gmres2,its,numR1solves,numR2solves] = GMRESAltSchwarz(tol,nmax,m,c,d,alpha,beta,gamma);
numR1solves
numR2solves
errorvgmres = norm(v_gmres - u_exact,inf)/norm(u_exact,inf)

subplot(3,4,10)
hold on
surf(X,Y,V_gmres1,'EdgeColor','none')
surf(X,Y,V_gmres2,'EdgeColor','none')
axis tight
colorbar
title('Case (v): GMRES w/ Alt. Schwarz')
hold off

%%

%%
% case (vi)
m=6; c=2/5; d=9/5; alpha=2; beta=7/2; gamma=2;
tol = 1e-9;
nmax = 100;

[u_exact,v_das,V_das1,V_das2,numR1solves,numR2solves,u1,u2] = DiscAltSchwarz(tol,nmax,m,c,d,alpha,beta,gamma);
numR1solves
numR2solves
errorvidas = norm(v_das - u_exact,inf)/norm(u_exact,inf)

x=h*(1:size(V_das1,1));
y=h*(1:size(V_das1,2));
[X,Y] = meshgrid(x,y);
X= X'; Y = Y';

subplot(3,4,11)
hold on
surf(X,Y,V_das1,'EdgeColor','none')
surf(X,Y,V_das2,'EdgeColor','none')
axis tight
colorbar
title('Case (vi): Disc. Alternating Schwarz')
hold off

[v_gmres,V_gmres1,V_gmres2,its,numR1solves,numR2solves] = GMRESAltSchwarz(tol,nmax,m,c,d,alpha,beta,gamma);
numR1solves
numR2solves
errorvigmres = norm(v_gmres - u_exact,inf)/norm(u_exact,inf)

subplot(3,4,12)
hold on
surf(X,Y,V_gmres1,'EdgeColor','none')
surf(X,Y,V_gmres2,'EdgeColor','none')
axis tight
colorbar
title('Case (vi): GMRES w/ Alt. Schwarz')
hold off

%%








