% Math 226B - Final Project
% Part 3b
% For each of the following two domain decomposition methods, write a 
% Matlab program that solves the discretized problem Av = b:
% (1) Multiplicative Schwarz method;
% The program for each of these methods should allow as inputs an arbitrary 
% initial guess v(0) for the solution of Av = b, a ?small? convergence 
% tolerance tol for the domain decomposition method, and a large integer 
% nmax to be used as a safeguard against an excessive number of iterations 
% of the domain decomposition method.


function [v_ms,V_ms] = MultiplicativeSchwarz(tol,nmax,m,c,d,alpha,beta,gamma)

[A,b,A1,A2,I1,I2,Num_Mat,tot_R1_pts,tot_R2_pts,tot_pts] = AbData(m,c,d,alpha,beta,gamma);

n = length(b);
h = 1/(2^m);

v0 = zeros(length(b));
v_whole = v0;

for n = 1:nmax
    
    if ((norm(A\b - v_whole)/norm(A\b - v0)) < tol)
        break;
    end
    
    w1 = b - A*v_whole;
    % reshape into a matrix
    W = zeros(size(Num_Mat,1),size(Num_Mat,2));
    for k = 1:tot_R1_pts
        [data_x,data_y] = find(Num_Mat == k);
        W(data_x,data_y) = w1(k);
    end
    solveA1 = fft2DPoisson(m,W(2:2^m,2:3*(2^m)),0,1,0,3);

    % reshape back into a vector
    solveA1_vec = zeros(n,1);
    for k = 1:tot_R1_pts
        [data_x,data_y] = find(Num_Mat == k);
        solveA1_vec(k) = solveA1(data_x-1,data_y-1);
    end
    v_half = v_whole + (I1.')*solveA1_vec; % requires solve with A1
    
    w2 = I2*(b - A*v_half);
    % reshape into a matrix
    R1R2_pts = tot_R1_pts + tot_R2_pts - tot_pts;
    for k = tot_R1_pts - R1R2_pts + 1:tot_pts
        [data_x,data_y] = find(Num_Mat == k);
        W(data_x,data_y) = w2(k);
    end
    solveA2 = fft2DPoisson(m,W((2^m)*c+1:end-1,(2^m)*d+1:(2^m)*d + 2^m -1)c,4,d,1);

    % reshape back into a vector
    solveA2_vec = zeros(n,1);
    for k = tot_R1_pts - R1R2_pts + 1:tot_pts
        [data_x,data_y] = find(Num_Mat == k);
        solveA2_vec(k) = solveA2(data_x,data_y);
    end

    v_whole = v_half + (I2.')*solveA2_vec; % reguires solve with A2
    
end

v_ms = v_whole; % this is a vector
n = length(b);
V_ms = zeros(size(Num_Mat,1),size(Num_Mat,2));

for i = 1:n
    [data_x,data_y] = find(Num_Mat == i);
    V_ms(data_x,data_y) = v_ms(i); 
end
    
    








% find relative error
% [X,Y]=meshgrid(0:h:c+4,0:h:3);
% X = X'; Y = Y';
% u = Y.^(alpha).*sin(beta*pi*X).*cos(gamma*pi*Y);
% 
% [max_V, max_ind] = max(V_ms(:));
% [row_ind,col_ind] = ind2sub(size(V_ms),max_ind);
% rel_error_ms = abs(max_V - u(row_ind,col_ind))/abs(u(row_ind,col_ind))



end



