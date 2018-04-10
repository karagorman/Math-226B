% Math 226B - Homework #3
% Problem 1c
% To test your Matlab program, use test cases that have solutions of the 
% form v(x, y) = y^? sin(??x) cos(??y), (2)
% where ? ? 0 and ?, ? > 0 are parameters. Determine the functions 
% f, b0, b1, c0, and c1 so that the function (2) is indeed the solution 
% of the above Poisson?s equation.

function [abs_error] = fft2DPoissonSolver(m,alpha,beta,gamma)

format long e

h = 1/(m+1);

x = [h:h:1-h];
y = [h:h:1-h];

[X,Y] = meshgrid(x,y);

% construct f matrix
f = (beta^2)*(pi^2).*(Y.^alpha).*sin(beta*pi.*X).*cos(gamma*pi.*Y)...
          + (pi^2)*(gamma^2).*(Y.^alpha).*sin(beta*pi.*X).*cos(gamma*pi.*Y)...
          + 2*pi*alpha*gamma.*(Y.^(alpha-1)).*sin(gamma*pi.*Y).*sin(beta*pi.*X)...
          - (alpha - 1)*alpha.*(Y.^(alpha - 2)).*cos(gamma*pi.*Y).*sin(beta*pi.*X);
        

% construct boundary conditions vectors
    if alpha == 0
        b0 = sin(beta*pi.*x);
    else
        b0 = zeros(1,m);
    end

b1 = sin(beta*pi.*x).*cos(gamma*pi);
c0 = zeros(1,m);
c1 = sin(beta*pi).*(y.^alpha).*cos((gamma*pi).*y);

% call Poisson solver with boundary conditions, and f as inputs
V = fft2DPoisson(m,b0,b1,c0,c1,f);

% construct exact solution matrix v_exact
v_exact = (Y.^alpha).*sin((beta*pi).*X).*cos((gamma*pi).*Y);

% compute absolute error
[max_V, max_ind] = max(V(:));
[row_ind, col_ind] = ind2sub(size(V),max_ind);
abs_error = abs(max_V - v_exact(row_ind,col_ind))


hold on

subplot(1,3,1)
mesh(X,Y,f);
xlabel('x');
ylabel('y');
zlabel('f(x,y)');
title('f(x,y)');


subplot(1,3,2);
mesh(X,Y,V);
xlabel('x');
ylabel('y');
zlabel('V(x,y)');
title('Approximate solution of V(x,y)')

subplot(1,3,3)
mesh(X,Y,v_exact);
xlabel('x');
ylabel('y');
zlabel('V(x,y)');
title('Exact V(x,y)');

end

