% function for MS

function [v_ms,V_ms,u,NumSolvesR1,NumSolvesR2] = MS(tol,nmax,m,c,d,...
    alpha,beta,gamma)

format long e

[A,b,A1,A2,I1,I2,Num_Mat,tot_R1_pts,tot_R2_pts,tot_pts,R1_xind,R1_yind,...
    R2_xind,R2_yind] = AbData(m,c,d,alpha,beta,gamma);

n = length(b);
h = 1/(2^m);

v0 = zeros(length(b),1);
v_whole = v0;
NumSolvesR1 = 0;
NumSolvesR2 = 0;

for nn = 1:nmax
    
    if ((norm(A\b - v_whole)/norm(A\b - v0)) < tol)
        nn
        break;
    end
    
    s = I1*(b - A*v_whole);
    %kmin = min(min(Num_Mat(R1_xind,R1_yind)));
    %kmax = max(max(Num_Mat(R1_xind,R1_yind)));
    kmin = 1;
    kmax = tot_R1_pts;
    
    w = RiSolve(R1_xind,R1_yind,kmin,kmax,s,0,1,0,3,Num_Mat,m);
    NumSolvesR1 = NumSolvesR1 + 1;
    
    v_half = v_whole + (I1.')*w;
    
    s = I2*(b - A*v_half);
    %kmin = min(min(Num_Mat(R2_xind,R2_yind)));
    %kmax = max(max(Num_Mat(R2_xind,R2_yind)));
    kmin = tot_pts - tot_R2_pts + 1;
    kmax = tot_pts;
    
    w = RiSolve(R2_xind,R2_yind,kmin,kmax,s,c,4,d,1,Num_Mat,m);
    NumSolvesR2 = NumSolvesR2 + 1;
    
    v_whole = v_half + (I2.')*w;
end

v_ms = v_whole; % this is a vector
n = length(b);
V_ms = zeros(size(Num_Mat,1),size(Num_Mat,2));

for i = 1:n
    [data_x,data_y] = find(Num_Mat == i);
    V_ms(data_x,data_y) = v_ms(i); 
end
    
% rel_res = norm(A*v_ms - b)/norm(b)
% 
x=0:h:c+4;
y=0:h:3;
[X,Y] = meshgrid(x,y);
X= X'; Y = Y';
% subplot(1,2,1)
% surf(X,Y,V_ms)
% axis tight
% colorbar
% title('Approx Solution with Multiplicative Schwarz')
% 
u = zeros(size(Num_Mat));
for i=1:n
   [data_x,data_y] = find(Num_Mat == i);
   u(data_x,data_y) = Y(data_x,data_y).^(alpha).*...
       sin(beta*pi*X(data_x,data_y)).*cos(gamma*pi*Y(data_x,data_y));
end
%    
% subplot(1,2,2)
% surf(X,Y,u)
% axis tight
% colorbar
% title('Exact Solution')

end