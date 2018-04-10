% Math 226B - Final Project
% Part 4a
% Create enumeration matrix for R1 and R2
% Inputs: m, a=c, b=d
% Output: 2 matrices. 1 corresponding to the enumeration of R1, and 1
% corresponding to the enumeration of R2
%clear all


function [EnumR1,EnumR2,justR2,tot_R1_pts,tot_Sigma1_pts,tot_Gamma1_pts,...
    totalR1pts,tot_R2_pts,tot_Sigma2_pts,tot_Gamma2_pts,totalR2pts,...
    R1_xind,R1_yind,R2_xind,R2_yind] = P4EnumeratedMats(m,c,d)

h = 1/(2^m);

x1 = 0:h:1;
y1 = 0:h:3;

x2 = c:h:c+4;
y2 = d:h:d+1;

nx1 = length(x1); 
ny1 = length(y1); 

nx2 = length(x2); 
ny2 = length(y2); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Build R1 Enumeration matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

EnumR1 = zeros(length(floor(0:h:(c+h+4))),ny1);

% R1 interior points
xLen_1 = 2^m - 1;
yLen_1 = length(h:h:3-h); 

vec_1 = 1:xLen_1*yLen_1;
mat_1 = reshape(vec_1,[xLen_1,yLen_1]);
EnumR1(2:xLen_1+1,2:yLen_1+1) = mat_1;

tot_R1_pts = length(vec_1);

% R1 boundary points = Sigma1
% bottom part of Simga1
xLen_Sigma1_1 = 2^m + 1;
yLen_Sigma1_1 = 1;

vec_Sigma1_1 = vec_1(end) + 1:vec_1(end) + xLen_Sigma1_1*yLen_Sigma1_1;
EnumR1(1:xLen_Sigma1_1,1) = vec_Sigma1_1.';

tot_Sigma1_pts = length(vec_Sigma1_1);

% left part of Sigma1
xLen_Sigma1_2 = 1;
yLen_Sigma1_2 = ny1-2;

vec_Sigma1_2 = vec_Sigma1_1(end) + 1:vec_Sigma1_1(end) + ...
    xLen_Sigma1_2*yLen_Sigma1_2;
EnumR1(1,2:yLen_Sigma1_2+1) = vec_Sigma1_2;

tot_Sigma1_pts = tot_Sigma1_pts + length(vec_Sigma1_2);

% top part of Sigma1
xLen_Sigma1_3 = 2^m + 1;
yLen_Sigma1_3 = 1;

vec_Sigma1_3 = vec_Sigma1_2(end) + 1:vec_Sigma1_2(end) + ...
    xLen_Sigma1_3*yLen_Sigma1_3;
EnumR1(1:xLen_Sigma1_3,end) = vec_Sigma1_3.';

tot_Sigma1_pts = tot_Sigma1_pts + length(vec_Sigma1_3);

% bottom right part of Sigma1
xLen_Sigma1_4 = 1;
yLen_Sigma1_4 = length(floor(h:h:d+h)) - 1;

vec_Sigma1_4 = vec_Sigma1_3(end) + 1:vec_Sigma1_3(end) + ...
    xLen_Sigma1_4*yLen_Sigma1_4;
EnumR1(2^m + 1,2:yLen_Sigma1_4+1) = vec_Sigma1_4;

tot_Sigma1_pts = tot_Sigma1_pts + length(vec_Sigma1_4);

% top right part of Sigma1
xLen_Sigma1_5 = 1;
yLen_SIgma1_5 = length(floor(3-h:-h:(d+1-h))) - 1;

vec_Sigma1_5 = vec_Sigma1_4(end) + 1:vec_Sigma1_4(end) + ...
    xLen_Sigma1_5*yLen_SIgma1_5;
EnumR1(2^m + 1,yLen_Sigma1_4+2^m +2:yLen_Sigma1_4+2^m + 1 + ...
    yLen_SIgma1_5) = vec_Sigma1_5;

tot_Sigma1_pts = tot_Sigma1_pts + length(vec_Sigma1_5);

% gamma1
xLen_Gamma1 = 1;
yLen_Gamma1 = 2^m;

vec_Gamma1 = vec_Sigma1_5(end) + 1:vec_Sigma1_5(end) + ...
    xLen_Gamma1*yLen_Gamma1;
EnumR1(2^m + 1,yLen_Sigma1_4 + 2:yLen_Sigma1_4 + 1 + ...
    yLen_Gamma1) = vec_Gamma1;

tot_Gamma1_pts = length(vec_Gamma1);

totalR1pts = tot_R1_pts + tot_Sigma1_pts + tot_Gamma1_pts;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Build R2 Enumeration matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

EnumR2 = zeros(nx2,ny2);
% R2 interior points
xLen_2 = length(h:h:4-h);
yLen_2 = 2^m - 1;

vec_2 = 1:xLen_2*yLen_2;
mat_2 = reshape(vec_2,[xLen_2,yLen_2]);
EnumR2(2:end-1,2:end-1) = mat_2;

tot_R2_pts = length(vec_2);

% bottom R2 boundary points = Sigma2
xLen_Sigma2_1 = nx2-length(c:h:1);
yLen_Sigma2_1 = 1;

vec_Sigma2_1 = vec_2(end) + 1:vec_2(end) + xLen_Sigma2_1*yLen_Sigma2_1;
EnumR2(length(c:h:1) + 1:end,1) = vec_Sigma2_1.';

tot_Sigma2_pts = length(vec_Sigma2_1);

% top part of Sigma2
xLen_Sigma2_2 = nx2-length(c:h:1);
yLen_Sigma2_2 = 1;

vec_Sigma2_2 = vec_Sigma2_1(end) + 1:vec_Sigma2_1(end) + ...
    xLen_Sigma2_2*yLen_Sigma2_2;
EnumR2(length(c:h:1) + 1:end,end) = vec_Sigma2_2.';

tot_Sigma2_pts = tot_Sigma2_pts + length(vec_Sigma2_2);

% right part of Sigma2
xLen_Sigma2_3 = 1;
yLen_Sigma2_3 = 2^m - 1;

vec_Sigma2_3 = vec_Sigma2_2(end) + 1:vec_Sigma2_2(end) + ...
    xLen_Sigma2_3*yLen_Sigma2_3;
EnumR2(end,2:end-1) = vec_Sigma2_3;

tot_Sigma2_pts = tot_Sigma2_pts + length(vec_Sigma2_3);

% bottom part of Gamma2
xLen_Gamma2_1 = length(c:h:1);
yLen_Gamma2_1 = 1;

vec_Gamma2_1 = vec_Sigma2_3(end) + 1:vec_Sigma2_3(end) + ...
    xLen_Gamma2_1*yLen_Gamma2_1;
EnumR2(1:xLen_Gamma2_1,1) = vec_Gamma2_1.';

tot_Gamma2_pts = length(vec_Gamma2_1);

% left part of Gamma2
xLen_Gamma2_2 = 1;
yLen_Gamma2_2 = 2^m - 1;

vec_Gamma2_2 = vec_Gamma2_1(end) + 1:vec_Gamma2_1(end) + ...
    xLen_Gamma2_2*yLen_Gamma2_2;
EnumR2(1,2:end-1) = vec_Gamma2_2;

tot_Gamma2_pts = tot_Gamma2_pts + length(vec_Gamma2_2);

% top part of Gamma2
xLen_Gamma2_3 = length(c:h:1);
yLen_Gamma2_3 = 1;

vec_Gamma2_3 = vec_Gamma2_2(end) + 1:vec_Gamma2_2(end) + ...
    xLen_Gamma2_3*yLen_Gamma2_3;
EnumR2(1:xLen_Gamma2_3,end) = vec_Gamma2_3.';

tot_Gamma2_pts = tot_Gamma2_pts + length(vec_Gamma2_3);
totalR2pts = tot_R2_pts + tot_Sigma2_pts + tot_Gamma2_pts;

justR2 = EnumR2;
% Pad with zeros to make EnumR2 the same size as EnumR1
EnumR2 = [zeros(2^m + 1 - length(c:h:1),size(EnumR2,2)); EnumR2];
EnumR2 = [zeros(size(EnumR2,1),floor(d/h)) , EnumR2, zeros(size(EnumR2,1),...
    3*(2^m) - (ceil((d+1)/h))+1)];

% Pad EnumR1 and EnumR2 with 0 on all 4 sides to represent going out of
% the grid
EnumR1 = [zeros(1,size(EnumR1,2)); EnumR1; zeros(1,size(EnumR1,2))];
EnumR1 = [zeros(size(EnumR1,1),1), EnumR1, zeros(size(EnumR1,1),1)];
EnumR2 = [zeros(1,size(EnumR2,2)); EnumR2; zeros(1,size(EnumR2,2))];
EnumR2 = [zeros(size(EnumR2,1),1), EnumR2, zeros(size(EnumR2,1),1)];

% probz not right. actually maybe right. check before using
R1_xind = 1 + (1:2^m + 1);
R1_yind = 1 + (1:ny1);
R2_xind = 1 + (2^m + 1 - length(c:h:1)+1:size(EnumR2,1)-2);
R2_yind = 1 + (floor(d/h)+1:floor(d/h)+2^m + 1);
end


