% Math 226B - Homework #5
% Problem 5
% Write a Matlab program that implements the Hermitian Lanczos process for 
% matrices A = AH ? Cn×n and starting vectors r ? Cn, r ?= 0, as presented i
% n class. Use an input parameter kmax to limit the maximum number of 
% Lanczos steps. Your program should only store the last two Lanczos 
% vectors, vk and vk?1, so that you can run the program for any number of 
% steps. The output of your routine should be the tridiagonal matrix 
% Tk ? Rk×k, where k is the iteration index at termination of your 
% algorithm.

function [Tk] = HermitianLanczos(kmax,m,fileNum)

format long e

if (fileNum == 1)
    load('HW5_P5a.mat')
elseif (fileNum == 2)
    load('HW5_P5b.mat')
end

A = make_3d_laplacian(m);
n = m^3;

tol = 1e-5;
alpha = zeros(kmax,1);
beta = zeros(kmax+1,1);

beta(1) = norm(r);
vk = r./beta(1);

for k=1:kmax
    
    q = A*vk;
    
    if (k > 1)
        q = q - beta(k)*v1;
    end
    
    alpha(k) = vk'*q;
    q = q - alpha(k)*vk;
    beta(k+1) = norm(q);
    
    if (beta(k+1) == 0)
        break
    end
    
    v1 = vk; %v1 = v_{k-1}
    vk=q./beta(k+1); %vk = v_{k+1}
end

% use alpha and beta vectors to construct Tk
Tk = diag(alpha,0) + diag(beta(2:kmax),-1) + diag(beta(2:kmax),1);

% now, use eig to compute the eigenvalues of Tk
lambda_approx = eig(Tk);
lambda_approx = uniquetol(lambda_approx,tol);

% compute exact eigenvalues to compare
i = [1:1:m]; j = [1:1:m]; l = [1:1:m];
[I,J,L] = meshgrid(i,j,l);
lambda_exact = 2*(3 - cos(I*pi/(m+1)) - cos(J*pi/(m+1)) - cos(L*pi/(m+1)));

lambda_exactVec = reshape(lambda_exact,[n,1]);
lambda_exactVec = uniquetol(lambda_exactVec,tol);

% for big case, find 10 smallest and 10 largest approximate and exact eigenvalues
if (fileNum == 2)
    lambda_approx = sort(lambda_approx);
    lambda_exactVec = sort(lambda_exactVec);
    
    smallest_approx = lambda_approx(1:10)
    smallest_exact = lambda_exactVec(1:10)
    
    biggest_approx = lambda_approx(end-9:end)
    biggest_exact = lambda_exactVec(end-9:end)
end
end
    
    
    
    
    
    
    