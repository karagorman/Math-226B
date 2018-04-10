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

[A,b,A1,A2,I1,I2,Num_Mat,tot_R1_pts,tot_R2_pts,tot_pts,R1_xind,R1_yind,R2_xind,R2_yind] = AbData(m,c,d,alpha,beta,gamma);

n = length(b);
h = 1/(2^m);

v0 = zeros(length(b),1);
v_whole = v0;

for n = 1:nmax
    
    if ((norm(A\b - v_whole)/norm(A\b - v0)) < tol)
        n
        break;
    end
    
    s = I1*(b - A*v_whole);
    % reshape into a matrix
    W = zeros(size(Num_Mat));
    for k = 1:tot_R1_pts
        [data_x,data_y] = find(Num_Mat == k);
        W(data_x,data_y) = s(k);
    end
    solveA1 = fft2DPoisson(m,W(2:2^m,2:3*(2^m)),0,1,0,3);
    %W(R1_xind,R1_yind) = fft2DPoisson(m,W(R1_xind,R1_yind),0,1,0,3);
    
    % pad with zeros to match Num_Mat
    solveA1 = [zeros(1,size(solveA1,2)); solveA1];
    solveA1 = [zeros(size(solveA1,1),1), solveA1];

    % reshape back into a vector
    w = zeros(length(s),1);
    for k = 1:tot_R1_pts
        [data_x,data_y] = find(Num_Mat == k);
        %w(k) = W(data_x,data_y);
        w(k) = solveA1(data_x,data_y);
    end
    
    v_half = v_whole + (I1.')*w; % requires solve with A1
    
    s = I2*(b - A*v_half);
    % reshape into a matrix
    R1R2_pts = tot_R1_pts + tot_R2_pts - tot_pts;
    kmin = min(min(Num_Mat(R2_xind,R2_yind)))
    kmax = max(max(Num_Mat(R2_xind,R2_yind)))
    %for k = tot_R1_pts - R1R2_pts + 1:tot_pts
    for k = 1:(kmax - kmin + 1)
        %[data_x,data_y] = find(Num_Mat == k);
        %W(data_x,data_y) = s(k - (tot_R1_pts - R1R2_pts + 1)+1);
        [data_x,data_y] = find(Num_Mat == kmin + k - 1);
        W(data_x,data_y) = s(k);
    end
    %solveA2 = fft2DPoisson(m,W((2^m)*c+1:end-1,(2^m)*d+1:(2^m)*d + 2^m -1),c,4,d,1);
    %solveA2 = fft2DPoisson(m,W(R2_xind,R2_yind),c,4,d,1);
  
    W(R2_xind,R2_yind) = fft2DPoisson(m,W(R2_xind,R2_yind),c,4,d,1);
    
    %solveA2 = [zeros((2^m)*c+1,size(solveA2,2)); solveA2];
    %solveA2 = [zeros(size(solveA2,1),(2^m)*d+1), solveA2]; 
    
    % reshape back into a vector
    w = zeros(length(s),1);
    for k = tot_R1_pts - R1R2_pts + 1:tot_pts
        [data_x,data_y] = find(Num_Mat == k);
        w(k-(tot_R1_pts - R1R2_pts + 1)+1) = W(data_x,data_y);
        %w(k-(tot_R1_pts - R1R2_pts + 1)+1) = solveA2(data_x,data_y);
    end

    v_whole = v_half + (I2.')*w; % reguires solve with A2
    
end

v_ms = v_whole; % this is a vector
n = length(b);
V_ms = zeros(size(Num_Mat,1),size(Num_Mat,2));

for i = 1:n
    [data_x,data_y] = find(Num_Mat == i);
    V_ms(data_x,data_y) = v_ms(i); 
end
    
rel_res = norm(A*v_ms - b)/norm(b)

x=0:h:c+4;
y=0:h:3;
[X,Y] = meshgrid(x,y);
X= X'; Y = Y';
subplot(1,2,1)
surf(X,Y,V_ms)
colorbar
title('Approx Solution with Multiplicative Schwarz')

u = zeros(size(Num_Mat));
for i=1:n
   [data_x,data_y] = find(Num_Mat == i);
   u(data_x,data_y) = Y(data_x,data_y).^(alpha).*sin(beta*pi*X(data_x,data_y)).*cos(gamma*pi*Y(data_x,data_y));
end
   
subplot(1,2,2)
surf(X,Y,u)
colorbar
title('Exact Solution')




% find relative error
% [X,Y]=meshgrid(0:h:c+4,0:h:3);
% X = X'; Y = Y';
% u = Y.^(alpha).*sin(beta*pi*X).*cos(gamma*pi*Y);
% 
% [max_V, max_ind] = max(V_ms(:));
% [row_ind,col_ind] = ind2sub(size(V_ms),max_ind);
% rel_error_ms = abs(max_V - u(row_ind,col_ind))/abs(u(row_ind,col_ind))



end



