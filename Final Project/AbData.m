% Math 226B - Final Project
% Part 3a
% Write Matlab programs that generate the ?data? A and b of the 
% discretized version, Av = b, of (1), as well as the submatrices A1 and 
% A2 of A corresponding to the subregions R1 and R2, respectively. Design 
% your program such that you can use any functions f : R  ? R and 
% g : ?R  ? R.


function [A,b,A1,A2,I1,I2,Num_Mat,tot_R1_pts,tot_R2_pts,tot_pts,R1_xind,...
    R1_yind,R2_xind,R2_yind] = AbData(m,c,d,alpha,beta,gamma)

if (0 < d && d < 2)
    [Num_Mat,tot_R1_pts,tot_R2_pts,tot_pts,R1_xind,R1_yind,R2_xind,...
        R2_yind] = NumberedMatd1(m,c,d);
elseif (d == 0)
    [Num_Mat,tot_R1_pts,tot_R2_pts,tot_pts,R1_xind,R1_yind,R2_xind,...
        R2_yind] = NumberedMatd0(m,c,d);
elseif (d== 2)
    [Num_Mat,tot_R1_pts,tot_R2_pts,tot_pts,R1_xind,R1_yind,R2_xind,...
        R2_yind] = NumberedMatd2(m,c,d);
end

h = 1/(2^m);

[X,Y]=meshgrid(0:h:c+4,0:h:3);
X = X'; Y = Y';

% construct f and g
% on the interior points.... from hw 3 problem 1
f = @(x,y) (beta^2)*(pi^2).*(y.^alpha).*sin(beta*pi.*x).*cos(gamma*pi.*y)...
          + (pi^2)*(gamma^2).*(y.^alpha).*sin(beta*pi.*x).*cos(gamma*pi.*y)...
          + 2*pi*alpha*gamma.*(y.^(alpha-1)).*sin(gamma*pi.*y).*sin(beta*pi.*x)...
          - (alpha - 1)*alpha.*(x.^(alpha - 2)).*cos(gamma*pi.*y).*sin(beta*pi.*x);

g = @(x,y) y.^(alpha).*sin(beta*pi*x).*cos(gamma*pi*y); % on the boundaries


% construct A and b
n = tot_pts;

A = 4*speye(n,n);
b = zeros(n,1);

for i=1:n
    [data_x,data_y] = find(Num_Mat == i);
    b(i) = h^2.*f(X(data_x,data_y),Y(data_x,data_y));
    
    % find the 4 closest points
    stencil = zeros(4,2);
    stencil(1,:) = [data_x, data_y + 1];
    stencil(2,:) = [data_x, data_y - 1];
    stencil(3,:) = [data_x + 1, data_y];
    stencil(4,:) = [data_x - 1, data_y];
    
    for j = 1:size(stencil,1)
        colInd = Num_Mat(stencil(j,1),stencil(j,2));
        if (colInd ~= 0) % then it is not a boundary
            A(i,colInd) = -1;
        elseif (colInd == 0) % then it is a boundary
            b(i) = b(i) + g(X(stencil(j,1),stencil(j,2)),...
                Y(stencil(j,1),stencil(j,2))); 
        end
    end
end



% now we need to get A1 and A2 from A
A1 = A(1:tot_R1_pts,1:tot_R1_pts);
A2 = A(tot_pts - tot_R2_pts + 1:end,tot_pts - tot_R2_pts + 1:end);

% construct I1 and I2
I1 = [speye(tot_R1_pts,tot_R1_pts), zeros(tot_R1_pts,tot_pts - tot_R1_pts)];
I2 = [zeros(tot_R2_pts,tot_pts - tot_R2_pts), speye(tot_R2_pts,tot_R2_pts)];

end






