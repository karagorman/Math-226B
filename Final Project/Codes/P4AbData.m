% Math 226B - Final Project
% Part 4a
% Generate the matrices A1, A2, B_gamma1, B_gamma2, IR2^gamma1, I_R1^gamma2
% and the vectors b1=h^2f1 - B_sigma1*g_sigma1, 
% b2=h^2f2 - B_sigma2*g_sigma2
% Inputs: EnumR1, EnumR2

function [EnumR1,EnumR2,justR2,tot_R1_pts,tot_Sigma1_pts,tot_Gamma1_pts,...
    totalR1pts,tot_R2_pts,tot_Sigma2_pts,tot_Gamma2_pts,totalR2pts,...
    R1_xind,R1_yind,R2_xind,R2_yind,A1,A2,B_Gamma1,B_Gamma2,b1,b2] ...
    = P4AbData(m,c,d,alpha,beta,gamma)

h = 1/(2^m);

[EnumR1,EnumR2,justR2,tot_R1_pts,tot_Sigma1_pts,tot_Gamma1_pts,totalR1pts,...
    tot_R2_pts,tot_Sigma2_pts,tot_Gamma2_pts,totalR2pts,R1_xind,R1_yind,...
    R2_xind,R2_yind] = P4EnumeratedMats(m,c,d);

% on the interior points.... from hw 3 problem 1
f = @(x,y) (beta^2)*(pi^2).*(y.^alpha).*sin(beta*pi.*x).*cos(gamma*pi.*y)...
          + (pi^2)*(gamma^2).*(y.^alpha).*sin(beta*pi.*x).*cos(gamma*pi.*y)...
          + 2*pi*alpha*gamma.*(y.^(alpha-1)).*sin(gamma*pi.*y).*sin(beta*pi.*x)...
          - (alpha - 1)*alpha.*(x.^(alpha - 2)).*cos(gamma*pi.*y).*sin(beta*pi.*x);

g = @(x,y) y.^(alpha).*sin(beta*pi*x).*cos(gamma*pi*y); % on the boundaries

X = -h:h:c+4+h;
Y = -h:h:3+h;

% generate matrices corresponding to R1
R1_mats = zeros(tot_R1_pts,totalR1pts);
f1 = zeros(tot_R1_pts,1);

for i=1:tot_R1_pts 
    [data_x,data_y] = find(EnumR1 == i);
    R1_mats(i,i) = 4;
    f1(i) = (h^2)*f(X(data_x),Y(data_y));
    
    % find the 4 closest points
    stencil = zeros(4,2);
    stencil(1,:) = [data_x, data_y + 1];
    stencil(2,:) = [data_x, data_y - 1];
    stencil(3,:) = [data_x + 1, data_y];
    stencil(4,:) = [data_x - 1, data_y];
    
    for j = 1:size(stencil,1)
        colInd = EnumR1(stencil(j,1),stencil(j,2));
        if (colInd ~= 0) % then it is not outside the boundary
            R1_mats(i,colInd) = -1;
        end
    end
end

% now we need to seperate into 3 matrices
A1 = R1_mats(1:tot_R1_pts,1:tot_R1_pts);
B_Sigma1 = R1_mats(1:tot_R1_pts,tot_R1_pts+1:tot_R1_pts+tot_Sigma1_pts);
B_Gamma1 = R1_mats(1:tot_R1_pts,tot_R1_pts+tot_Sigma1_pts+1:end);

g_Sigma1 = zeros(tot_Sigma1_pts,1);
for i = tot_R1_pts + 1:tot_R1_pts + tot_Sigma1_pts
    [data_x,data_y] = find(EnumR1 == i);
    g_Sigma1(i - tot_R1_pts) = g(X(data_x),Y(data_y));
end
b1 = f1 - B_Sigma1*g_Sigma1;

% generate matrices corresponding to R2
R2_mats = zeros(tot_R2_pts,totalR2pts);
f1 = zeros(tot_R2_pts,1);

for i=1:tot_R2_pts 
    [data_x,data_y] = find(EnumR2 == i);
    R2_mats(i,i) = 4;
    f2(i) = (h^2)*f(X(data_x),Y(data_y));
    
    % find the 4 closest points
    stencil = zeros(4,2);
    stencil(1,:) = [data_x, data_y + 1];
    stencil(2,:) = [data_x, data_y - 1];
    stencil(3,:) = [data_x + 1, data_y];
    stencil(4,:) = [data_x - 1, data_y];
    
    for j = 1:size(stencil,1)
        colInd = EnumR2(stencil(j,1),stencil(j,2));
        if (colInd ~= 0) % then it is not outside the boundary
            R2_mats(i,colInd) = -1;
        end
    end
end

% now we need to seperate into 3 matrices
A2 = R2_mats(1:tot_R2_pts,1:tot_R2_pts);
B_Sigma2 = R2_mats(1:tot_R2_pts,tot_R2_pts+1:tot_R2_pts+tot_Sigma2_pts);
B_Gamma2 = R2_mats(1:tot_R2_pts,tot_R2_pts+tot_Sigma2_pts+1:end);

g_Sigma2 = zeros(tot_Sigma2_pts,1);
for i = tot_R2_pts + 1:tot_R2_pts + tot_Sigma2_pts
    [data_x,data_y] = find(EnumR2 == i);
    g_Sigma2(i - tot_R2_pts) = g(X(data_x),Y(data_y));
end
b2 = f2.' - B_Sigma2*g_Sigma2;
end
    
    
    

