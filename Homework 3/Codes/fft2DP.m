function [abs_error] = fft2DP(m,alpha,beta,gamma)

format long e

h = 1/(m+1);
% n = m-2; % number of interior points
x = [h:h:1-h];
y = [h:h:1-h];

% construct f matrix and exact v matrix
% f = zeros(m,m);
v_exact = zeros(m,m);
for j = 1:m
    for k = 1:m
        f(j,k) = (beta^2)*(pi^2)*(y(k)^alpha)*sin(beta*pi*x(j))*cos(gamma*pi*y(k))...
                   + (pi^2)*(gamma^2)*(y(k)^alpha)*sin(beta*pi*x(j))*cos(gamma*pi*y(k))...
                   + 2*pi*alpha*gamma*(y(k)^(alpha-1))*sin(gamma*pi*y(k))*sin(beta*pi*x(j))...
                   - (alpha - 1)*alpha*(y(k)^(alpha - 2))*cos(gamma*pi*y(k))*sin(beta*pi*x(j));
        
        v_exact(j,k) = (y(k)^alpha)*sin(beta*pi*x(j))*cos(gamma*pi*y(k));
    end
end

% boundary conditions
    if alpha == 0
        b0 = sin(beta*pi*x);
    else
        b0 = zeros(1,m);
    end

b1 = sin(beta*pi*x).*cos(gamma*pi);
c0 = zeros(1,m);
c1 = (y.^alpha).*sin(beta*pi).*cos(gamma*pi*y);

f(1,:) = f(1,:) + 1/(h^2)*b0;
f(m,:) = f(m,:) + 1/(h^2)*b1;
f(:,1) = f(:,1) + 1/(h^2)*c0';
f(:,m) = f(:,m) + 1/(h^2)*c1';


f_bar = zeros(m,m);

for k = 1:m
    f_bar(:,k) = fftMult(f(:,k));
end

for j = 1:m
    f_bar(j,:) = fftMult(f_bar(j,:).').';
end


[X,Y]=meshgrid(1:m,1:m);
lambda = 2*(1 - cos(pi*h*X) + 1 - cos(pi*h*Y));
v_bar = (h^2)*f./lambda


for k = 1:m
    v(:,k) = fftMult(v_bar(:,k));
end

for j = 1:m
    v(j,:) = fftMult(v(j,:).').';
end




[max_v, max_ind] = max(v(:));
[row_ind, col_ind] = ind2sub(size(v),max_ind);
abs_error = abs(max_v - v_exact(row_ind,col_ind))


[X, Y] = meshgrid(x,y);

hold on

subplot(1,3,1)
mesh(X,Y,f)
xlabel('x')
ylabel('y')
zlabel('f')
title('f(x,y) for 2D Poisson Equation')

% x1 = h*(1:m);
% y1 = h*(1:m);
% [X1, Y1] = meshgrid(x1,y1);

subplot(1,3,2)
mesh(X,Y,v)
xlabel('x')
ylabel('y')
zlabel('V')
title('Approximation of v(x,y) for 2D Poisson Equation')

subplot(1,3,3)
mesh(X,Y,v_exact)
xlabel('x')
ylabel('y')
zlabel('v_exact')
title('Exact v(x,y) for 2D Poisson Equation')

    function w = fftMult(x)
        
        n = length(x);
        x_tilde = [0; x; zeros(n+1,1)];
        w_tilde = fft(x_tilde);
        w_hat = w_tilde(2:n+1);
        w = -sqrt(2*h)*imag(w_hat);
        
    end
end
