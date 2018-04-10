% Math 226B - Homework #5
% Problem 4
% Write a Matlab program that implements the Arnoldi process for matrices 
% A ? Cn×n and starting vectors r ? Cn, r ?= 0, as presented in class. U
% se an input parameter kmax to limit the maximum number of Arnoldi steps. 
% The output of your program should be the upper-Hessenberg matrix 
% Hk ? Ck×k and the matrix Vk ? Cn×k with the first k Arnoldi vectors, 
% where k is the iteration index at termination of your algorithm.

function [H,V] = Arnoldi(kmax)

format long e

load('HW5_P4.mat') % uploads r vector

m = 100;
n = m^2;
k = kmax;
beta = norm(r);
V = zeros(n,k+1);
V(:,1) = r/beta;
H_tilde = zeros(k+1,k);

gammaVec = [1, 10, 50, 100, 1000];

for g = 1:length(gammaVec)
    for j = 1:k
        %q = A*V(:,j);
        q = ApMult(V(:,j),gammaVec(g),m);

        for i = 1:j
            H_tilde(i,j) = V(:,i)'*q;
            q = q - H_tilde(i,j)*V(:,i);
        end

        H_tilde(j+1,j) = norm(q);

        if (H_tilde(j+1,j) == 0)
            k = j;
            break
        end

        V(:,j+1) = q/H_tilde(j+1,j);

    end

    % compute eigenpairs of Hk
    Hk = H_tilde(1:k,1:k);
    [z,lambdaMat] = eig(Hk);
    lambda_approx = diag(lambdaMat);

    % compute the residual norm of the approximate eigenpairs without 
    % computing the approximate eigenvectors
    % rho = zeros(k,1);
    j = [1:1:k];
    rho = H_tilde(k+1,k).*abs(z(k,j));
    rho = rho';
    gammaVec(g)
    error = (norm(rho))

    min_rho_val = min(rho)
    max_rho_val = max(rho)
    [row_min] = find(rho == min_rho_val);
    [row_max] = find(rho == max_rho_val);
    
    row_min
    min_lams = lambda_approx(row_min)
    row_max
    max_lams = lambda_approx(row_max)

    % plot that shows all k approximate eigenvalues of A
    figure(1)
    hold on
    subplot(2,3,g)
    plot(1:k,lambda_approx)
    xlabel('k')
    ylabel("Approximate Eigenvalue, \lambda_k, of A")
    title(strcat("\gamma = ", num2str(gammaVec(g))))
    
    figure(2)
    subplot(2,3,g)
    plot(real(lambda_approx),imag(lambda_approx),'x')
    xlabel("Real Part of \lambda_k")
    ylabel("Imaginary Part of \lambda_k")
    title(strcat("\gamma = ", num2str(gammaVec(g))))
end
end