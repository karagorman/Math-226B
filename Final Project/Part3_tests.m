% To Run Part 3
%

clear all



%%
% case(i): c=1/2, d=0, alpha=0, beta=1, gamma=1/2
tic;
m = 6; c = 1/2; d = 0; alpha = 0; beta = 1; gamma = 1/2; %m=6
h = 1/(2^m);
tol = 1e-6;
nmax = 100;
[v_ms,V_ms,u,NumSolvesR1,NumSolvesR2] = MS(tol,nmax,m,c,d,alpha,beta,gamma);
NumSolvesR1
NumSolvesR2

% plot ms
x=0:h:c+4;
y=0:h:3; 
[X,Y] = meshgrid(x,y);
X= X'; Y = Y';

errorims = max(max(abs(V_ms - u)))/max(max(abs(u)))

subplot(3,4,1) 
surf(X,Y,V_ms,'Edgecolor','none')
colorbar
title('Case (i): Multiplicative Schwarz')

[v_pcgs,V_pcgs,u,NumSolvesR1,NumSolvesR2] = PCGSchwarz(tol,nmax,m,c,d,...
    alpha,beta,gamma);
NumSolvesR1
NumSolvesR2
erroripcgs = max(max(abs(V_pcgs - u)))/max(max(abs(u)))

%plot pcgs
subplot(3,4,2)
surf(X,Y,V_ms,'Edgecolor','none')
colorbar
title('Case (i): PCG Schwarz')
toc;
%%




%%
% case(ii): c=1/2, d=0, alpha=2, beta=7/2, gamma=2 %alex
tic;
m = 7; c = 1/2; d = 0; alpha = 2; beta = 7/2; gamma = 2; %m=7 probz 
h = 1/(2^m);
tol = 1e-6;
nmax = 100;
[v_ms,V_ms,u,NumSolvesR1,NumSolvesR2] = MS(tol,nmax,m,c,d,alpha,beta,gamma);
NumSolvesR1
NumSolvesR2
% plot ms
x=0:h:c+4;
y=0:h:3; 
[X,Y] = meshgrid(x,y);
X= X'; Y = Y';

erroriims = max(max(abs(V_ms - u)))/max(max(abs(u)))

subplot(3,4,3) 
surf(X,Y,V_ms,'Edgecolor','none')
colorbar
title('Case (ii): Multiplicative Schwarz')

m=6;
[v_pcgs,V_pcgs,u,NumSolvesR1,NumSolvesR2] = PCGSchwarz(tol,nmax,m,c,d,alpha,beta,gamma);
NumSolvesR1
NumSolvesR2
erroriipcgs = max(max(abs(V_pcgs - u)))/max(max(abs(u)))

%plot pcgs
subplot(3,4,4)
surf(X,Y,V_ms,'Edgecolor','none')
colorbar
title('Case (ii): PCG Schwarz')
toc;
%%


%%
% case(iii): c=1/2, d=1, alpha=0, beta=1, gamma=1/2 %adam
tic;
m = 6; c = 1/2; d = 1; alpha = 0; beta = 1; gamma = 1/2;
h = 1/(2^m);
tol = 1e-6;
nmax = 100;
[v_ms,V_ms,u,NumSolvesR1,NumSolvesR2] = MS(tol,nmax,m,c,d,alpha,beta,gamma);
NumSolvesR1
NumSolvesR2
% plot ms
x=0:h:c+4;
y=0:h:3; 
[X,Y] = meshgrid(x,y);
X= X'; Y = Y';

erroriiims = max(max(abs(V_ms - u)))/max(max(abs(u)))

subplot(3,4,5) 
surf(X,Y,V_ms,'Edgecolor','none')
colorbar
title('Case (iii): Multiplicative Schwarz')

[v_pcgs,V_pcgs,u,NumSolvesR1,NumSolvesR2] = PCGSchwarz(tol,nmax,m,c,d,alpha,beta,gamma);
NumSolvesR1
NumSolvesR2
erroriiipcgs = max(max(abs(V_pcgs - u)))/max(max(abs(u)))

%plot pcgs
subplot(3,4,6)
surf(X,Y,V_ms,'Edgecolor','none')
colorbar
title('Case (iii): PCG Schwarz')
toc;
%%



%%
% case(iv): c=1/2, d=1, alpha=2, beta=7/2, gamma=2 %adam
tic;
m = 7; c = 1/2; d = 1; alpha = 2; beta = 7/2; gamma = 2;
h = 1/(2^m);
tol = 1e-6;
nmax = 100;
[v_ms,V_ms,u,NumSolvesR1,NumSolvesR2] = MS(tol,nmax,m,c,d,alpha,beta,gamma);
NumSolvesR1
NumSolvesR2
% plot ms
x=0:h:c+4;
y=0:h:3; 
[X,Y] = meshgrid(x,y);
X= X'; Y = Y';

errorivms = max(max(abs(V_ms - u)))/max(max(abs(u)))

subplot(3,4,7) 
surf(X,Y,V_ms,'Edgecolor','none')
colorbar
title('Case (iv): Multiplicative Schwarz')

m = 5;
[v_pcgs,V_pcgs,u,NumSolvesR1,NumSolvesR2] = PCGSchwarz(tol,nmax,m,c,d,alpha,beta,gamma);
NumSolvesR1
NumSolvesR2
errorivpcgs = max(max(abs(V_pcgs - u)))/max(max(abs(u)))

%plot pcgs
subplot(3,4,8)
surf(X,Y,V_ms,'Edgecolor','none')
colorbar
title('Case (iv): PCG Schwarz')
toc;
%%


%%
% case(v); c=1/2, d=2, alpha=0, beta=1, gamma=1/2 %me
tic;
m = 6; c = 1/2; d = 2; alpha = 0; beta = 1; gamma = 1/2;
h = 1/(2^m);
tol = 1e-6;
nmax = 100;
[v_ms,V_ms,u,NumSolvesR1,NumSolvesR2] = MS(tol,nmax,m,c,d,alpha,beta,gamma);
NumSolvesR1
NumSolvesR2
% plot ms
x=0:h:c+4;
y=0:h:3;
[X,Y] = meshgrid(x,y);
X= X'; Y = Y';

errorvms = max(max(abs(V_ms - u)))/max(max(abs(u)))

subplot(3,4,9) 
surf(X,Y,V_ms,'Edgecolor','none')
colorbar
title('Case (v): Multiplicative Schwarz')

nmax = 20;
[v_pcgs,V_pcgs,u,NumSolvesR1,NumSolvesR2] = PCGSchwarz(tol,nmax,m,c,d,alpha,beta,gamma);
NumSolvesR1
NumSolvesR2
errorvpcgs = max(max(abs(V_pcgs - u)))/max(max(abs(u)))

%plot pcgs
subplot(3,4,10)
surf(X,Y,V_ms,'Edgecolor','none')
colorbar
title('Case (v): PCG Schwarz')
toc;
%%




%%
% case(vi): c=1/2, d=2, alpha=2, beta=7/2, gamma=2 % me
tic;
m = 7; c = 1/2; d = 2; alpha = 2; beta = 7/2; gamma = 2; %m=7 probz, but not for cg
h = 1/(2^m);
tol = 1e-6;
nmax = 100;
[v_ms,V_ms,u,NumSolvesR1,NumSolvesR2] = MS(tol,nmax,m,c,d,alpha,beta,gamma);
NumSolvesR1
NumSolvesR2
% plot ms
x=0:h:c+4;
y=0:h:3;
[X,Y] = meshgrid(x,y);
X= X'; Y = Y';

errorvims = max(max(abs(V_ms - u)))/max(max(abs(u)))

subplot(3,4,11) 
surf(X,Y,V_ms,'Edgecolor','none')
colorbar
title('Case (vi): Multiplicative Schwarz')

m=6; nmax=25;
[v_pcgs,V_pcgs,u,NumSolvesR1,NumSolvesR2] = PCGSchwarz(tol,nmax,m,c,d,alpha,beta,gamma);
NumSolvesR1
NumSolvesR2
errorvipcgs = max(max(abs(V_pcgs - u)))/max(max(abs(u)))

%plot pcgs
subplot(3,4,12)
surf(X,Y,V_ms,'Edgecolor','none')
colorbar
title('Case (vi): PCG Schwarz')
toc;
%%





