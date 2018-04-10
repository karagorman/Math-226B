% Math 226B - Final Project
% Part 4b
% GMRES with Alternating Schwarz as Preconditioner


function [v_gmres,V_gmres1,V_gmres2,its,numR1solves,numR2solves] = ...
    GMRESAltSchwarz(tol,nmax,m,c,d,alpha,beta,gamma)

format long e

[EnumR1,EnumR2,justR2,tot_R1_pts,tot_Sigma1_pts,tot_Gamma1_pts,totalR1pts,...
    tot_R2_pts,tot_Sigma2_pts,tot_Gamma2_pts,totalR2pts,R1_xind,R1_yind,...
    R2_xind,R2_yind,A1,A2,B_Gamma1,B_Gamma2,b1,b2] = ...
    P4AbData(m,c,d,alpha,beta,gamma);

h = 1/(2^m);
b = [b1;b2];
n = length(b);
v0 = zeros(length(b),1);

v0p = v0;
R1Interior = EnumR1(3:2^m +1,3:end-2);
R2Interior = justR2(2:end-1,2:end-1);
t1 = RiSolve(1:size(R1Interior,1),1:size(R1Interior,2),1,tot_R1_pts,b1.',...
    0,1,0,3,R1Interior,m);
t2 = RiSolve(1:size(R2Interior,1),1:size(R2Interior,2),1,tot_R2_pts,b2,...
    c,4,d,1,R2Interior,m);
bp=[t1;t2];
numR1solves = 1;
numR2solves = 1;


[v_gmres,flag,relres,iter] = gmres(@(v) ApMult(v),bp,[],tol,nmax,[],[],v0p);
its = iter(2);
numR1solves = numR1solves + its;
numR2solves = numR2solves + its;

V_gmres1 = zeros(size(EnumR1));
V_gmres2 = zeros(size(EnumR2));

for i = 1:tot_R1_pts
    [data_x,data_y] = find(EnumR1 == i);
    %V_gmres1(data_x,data_y) = v1(i);
    V_gmres1(data_x,data_y) = v_gmres(i);
end

for i = 1:tot_R2_pts
    [data_x,data_y] = find(EnumR2 == i);
    V_gmres2(data_x,data_y) = v_gmres(i+tot_R1_pts);
end


% x=h*(1:size(EnumR1,1));
% y=h*(1:size(EnumR1,2));
% [X,Y] = meshgrid(x,y);
% X= X'; Y = Y';
% hold on
% surf(X,Y,V_gmres1,'EdgeColor','none')
% surf(X,Y,V_gmres2,'EdgeColor','none')
% axis tight
% colorbar

    % function to compute q=A'p
    function q = ApMult(p)
         p1 = p(1:size(A1,2),1);
         p2 = p(length(p1)+1:end,1);

         p2_bar = InterpolationR2Gamma1(m,EnumR1,EnumR2,p2,tot_R1_pts,...
             tot_Sigma1_pts,totalR1pts);
         R1Interior = EnumR1(3:2^m +1,3:end-2);
         t1 = RiSolve(1:size(R1Interior,1),1:size(R1Interior,2),1,...
             tot_R1_pts,B_Gamma1*p2_bar,0,1,0,3,R1Interior,m); 
         t1_bar = InterpolationR1Gamma2(m,EnumR1,EnumR2,t1,tot_R2_pts,...
             tot_Sigma2_pts,totalR2pts);
         R2Interior = justR2(2:end-1,2:end-1);
         t2 = RiSolve(1:size(R2Interior,1),1:size(R2Interior,2),1,...
             tot_R2_pts,B_Gamma2*t1_bar,c,4,d,1,R2Interior,m); 
         q = [p1+t1; p2 - t2];
    end
     
end 






