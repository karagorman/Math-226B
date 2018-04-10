% Math 226B - Homework #5
% Problem 1




function [u,relres,resvec] = Poisson2DGMRES(m,gamma)

format long e

n = m^2;

h = 1/(m+1);

x = h:h:1-h;
y = h:h:1-h;
[X,Y] = meshgrid(x,y);
X = X';
Y = Y';

% construct f matrix
f = (X.^3).*(Y.^2).*exp(2 - X - Y);

% construct boundary condition vectors based on g(x,y)
b0 = ones(1,m);
b1 = ones(1,m);
c0 = ones(m,1);
c1 = ones(m,1);

b = h^2.*f;

b(1,:) = b(1,:) + b0;
b(m,:) = b(m,:) + b1;
b(:,1) = b(:,1) + gamma*(h/2)*c0 + c0;
b(:,m) = b(:,m) - gamma*(h/2)*c1 + c1;

b = reshape(b,[m,m]);
bp = fft2DPoisson(m,b);
bp = reshape(bp,[n,1]);

tol = 1e-10;
maxit = n;

% run GMRES
[u,flag,relres,iter,resvec] = gmres(@(z) ApMult(z,gamma,m),bp,[],tol,maxit);
total_its = iter(2)
relres
resvec = resvec./resvec(1);
end





