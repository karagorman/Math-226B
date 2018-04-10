% Math 226B - Final Project
% Part 3a
% Function to created enumerated matrix if b=1 (but more specifically, if b
% does not equal 0 or 2)
% where we are renaming b to d because b is also the name of the right hand
% side of Av = b
% inputs: m, c, d
% outputs: Num_Mat, vec_5 (vec_5(end) will be the size of matrix A, i.e.,
% the number of points in the enumerated matrix)

function [Num_Mat,tot_R1_pts,tot_R2_pts,tot_pts,R1_xind,R1_yind,R2_xind,...
    R2_yind] = NumberedMatd1(m,c,d)

h = 1/(2^m);

x = h:h:c+4-h;
y = h:h:3-h;

nx = length(x);
ny = length(y);

Num_Mat = zeros(nx,ny);

% part of R1 below R2
xLen_1 = 2^m - 1;
yLen_1 = (2^m)*d-1;

vec_1 = 1:xLen_1*yLen_1;
mat_1 = reshape(vec_1,[xLen_1, yLen_1]);
Num_Mat(1:xLen_1,1:yLen_1) = mat_1;

R1_pts = length(vec_1);

% part of R1 to the left of R2
xLen_2 = (2^m)*c - 1;
yLen_2 = 2^m + 1;

vec_2 = vec_1(end)+1:vec_1(end) + xLen_2*yLen_2;
mat_2 = reshape(vec_2,[xLen_2,yLen_2]);
Num_Mat(1:xLen_2,yLen_1 + 1:yLen_1 + yLen_2) = mat_2;

R1_pts = R1_pts + length(vec_2);

% part of R1 above R2
xLen_3 = 2^m - 1;
yLen_3 = (2^m)*(2-d) - 1;

vec_3 = vec_2(end) + 1:vec_2(end) + xLen_3*yLen_3;
mat_3 = reshape(vec_3,[yLen_3,xLen_3]);
Num_Mat(1:xLen_3,yLen_1+ yLen_2 + 1:yLen_1 + yLen_2 + yLen_3) = mat_3;

R1_pts = R1_pts + length(vec_3);

% first part of gamma2
xLen_gamma2_1 = (2^m - 1) - ((2^m)*c - 1);
yLen_gamma2_1 = 1;

vec_gamma2_1 = vec_3(end) + 1:vec_3(end) + xLen_gamma2_1*yLen_gamma2_1;
Num_Mat(xLen_2 + 1:xLen_1,2^m*d) = vec_gamma2_1';

gamma2_pts = length(vec_gamma2_1);

% second part of gamma2
xLen_gamma2_2 = 1;
yLen_gamma2_2 = 2^m - 1;

vec_gamma2_2 = vec_gamma2_1(end) + 1:vec_gamma2_1(end) + ...
    xLen_gamma2_2*yLen_gamma2_2;
Num_Mat(xLen_2+1,yLen_1 + 2:yLen_1 +1 + yLen_gamma2_2) = vec_gamma2_2;

gamma2_pts = gamma2_pts + length(vec_gamma2_2);

% third part of gamma2
xLen_gamma2_3 = (2^m - 1) - ((2^m)*c - 1);
yLen_gamma2_3 = 1;

vec_gamma2_3 = vec_gamma2_2(end) + 1:vec_gamma2_2(end) + ...
    xLen_gamma2_3*yLen_gamma2_3;
Num_Mat(xLen_2 + 1:xLen_1,(2^m)*(d+1)) = vec_gamma2_3';

gamma2_pts = gamma2_pts + length(vec_gamma2_3);

% intersecting part of R1 and R2
xLen_4 = (2^m)*(1-c) - 1;
yLen_4 = 2^m - 1;

vec_4 = vec_gamma2_3(end) + 1:vec_gamma2_3(end) +xLen_4*yLen_4;
mat_4 = reshape(vec_4,[xLen_4,yLen_4]);
Num_Mat((2^m)*c+1:2^m-1,(2^m*d + 1):(2^m)*(d+1) - 1) = mat_4;

R1R2_pts = length(vec_4);

% gamma1
xLen_gamma1 = 1;
yLen_gamma1 = 2^m - 1;

vec_gamma1 = vec_4(end) + 1:vec_4(end) + xLen_gamma1*yLen_gamma1;
Num_Mat(2^m,(2^m*d + 1):(2^m)*(d+1) - 1) = vec_gamma1;

gamma1_pts = length(vec_gamma1);

% part of R2 sticking out
xLen_5 = (2^m)*(c+3) - 1;
yLen_5 = 2^m - 1;

vec_5 = vec_gamma1(end) + 1:vec_gamma1(end) + xLen_5*yLen_5;
mat_5 = reshape(vec_5,[xLen_5,yLen_5]);
Num_Mat(2^m + 1:end,(2^m)*d + 1:(2^m)*(d+1)-1) = mat_5;

R2_pts = length(vec_5);

% put zeros on boundary of Num_Mat
Num_Mat = [zeros(size(Num_Mat,1),1), Num_Mat , zeros(size(Num_Mat,1),1)];
Num_Mat = [zeros(1,size(Num_Mat,2)); Num_Mat ; zeros(1,size(Num_Mat,2))];

% total number of points
tot_R1_pts = R1_pts + gamma2_pts + R1R2_pts;
tot_R2_pts = R1R2_pts + gamma1_pts + R2_pts;
tot_pts = tot_R1_pts + tot_R2_pts - R1R2_pts;

% make vectors of indices of R1 and R2 points
R1_xind = 1 + (1:2^m - 1);
R1_yind = 1 + (1:ny);
R2_xind = 1 + ((2^m)*c+1:nx);
R2_yind = 1 + ((2^m)*d+1:(2^m)*(d+1)-1);

end