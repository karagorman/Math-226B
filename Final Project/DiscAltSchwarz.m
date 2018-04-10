% Math 226B - Final Project
% Part 4b
% Discretized Alternating Schwarz Method

function [u_exact,v_das,V_das1,V_das2,numR1solves,numR2solves,u1,u2] = ...
    DiscAltSchwarz(tol,nmax,m,c,d,alpha,beta,gamma)

format long e

[EnumR1,EnumR2,justR2,tot_R1_pts,tot_Sigma1_pts,tot_Gamma1_pts,totalR1pts,...
    tot_R2_pts,tot_Sigma2_pts,tot_Gamma2_pts,totalR2pts,R1_xind,R1_yind,...
    R2_xind,R2_yind,A1,A2,B_Gamma1,B_Gamma2,b1,b2] = P4AbData(m,c,d,...
    alpha,beta,gamma);

h = 1/(2^m);
b = [b1;b2];
n = length(b);
v2 = zeros(size(A2,1),1);


for n = 1:nmax
    v_Gamma1 = InterpolationR2Gamma1(m,EnumR1,EnumR2,v2,tot_R1_pts,...
        tot_Sigma1_pts,totalR1pts);
    b1_bar = b1 - B_Gamma1*v_Gamma1; 
    kmin = 1;
    kmax = tot_R1_pts;
    R1Interior = EnumR1(3:2^m +1,3:end-2);
    v1 = RiSolve(1:size(R1Interior,1),1:size(R1Interior,2),kmin,kmax,...
        b1_bar,0,1,0,3,R1Interior,m);
    
    v_Gamma2 = InterpolationR1Gamma2(m,EnumR1,EnumR2,v1,tot_R2_pts,...
        tot_Sigma2_pts,totalR2pts);
    b2_bar = b2 - B_Gamma2*v_Gamma2;
    kmin = 1;
    kmax = tot_R2_pts;
    R2Interior = justR2(2:end-1,2:end-1);
    v2 = RiSolve(1:size(R2Interior,1),1:size(R2Interior,2),kmin,kmax,...
        b2_bar,c,4,d,1,R2Interior,m);
    
    vn = [v1; v2];
    Avn = [A1*v1 + B_Gamma1*v_Gamma1; B_Gamma2*v_Gamma2 + A2*v2];
    if (norm(b - Avn) < tol)
        numR1solves = n;
        numR2solves = n;
        break
    end
end
v_das = vn;
V_das1 = zeros(size(EnumR1));
V_das2 = zeros(size(EnumR2));

for i = 1:length(v1)
    [data_x,data_y] = find(EnumR1 == i);
    V_das1(data_x,data_y) = v1(i); 
end

for i = 1:length(v2)
    [data_x,data_y] = find(EnumR2 == i);
    V_das2(data_x,data_y) = v2(i); 
end


x=h*(1:size(V_das1,1));
y=h*(1:size(V_das1,2));
[X,Y] = meshgrid(x,y);
X= X'; Y = Y';
u = (Y.^alpha).*sin(beta*pi*X).*cos(gamma*pi*Y);

u1 = zeros(length(v1),1);
for i=1:length(v1)
   [data_x,data_y] = find(EnumR1 == i);
    u1(i) = u(data_x,data_y);
end

u2 = zeros(length(v2),1);
for i=1:length(v2)
   [data_x,data_y] = find(EnumR2 == i);
    u2(i) = u(data_x,data_y);
end


u_exact = [u1; u2];
% error = norm(v_das - u_exact,inf)/norm(u_exact,inf)



% hold on
% %subplot(1,2,1)
% surf(X,Y,V_das1,'EdgeColor','none')
% %surf(X,Y,V_das1)
% %axis tight
% %colorbar
% %title('Approx Solution of R1 with Disc Alt Schwarz')
% 
% %subplot(1,2,2)
% surf(X,Y,V_das2,'EdgeColor','none')
% %surf(X,Y,V_das2)
% axis tight
% colorbar
% %title('Approx Solution of R2 with Disc Alt Schwarz')


end






