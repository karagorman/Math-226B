% Math 226B - Homework #1
% Problem 2
% Formulate a linear algebra problem the solution of which is the PageRank 
% vector x of this internal network and compute x.

%n = 10;

%Q = [0,0,0,0,1/2,0,1/2,0,0,0;
%     0,0,0,0,0,1/3,1/3,0,1/3,0;
%     0,0,0,0,0,0,0,0,0,0;
%     0,0,0,0,0,0,1,0,0,0;
%     0,0,0,0,0,0,1/3,0,1/3,1/3;
%     0,0,0,0,0,0,0,0,0,0;
%     0,0,0,0,0,0,0,0,0,0;
%     0,0,1/4,1/4,0,0,1/4,0,1/4,0;
%     0,1/5,0,1/5,0,0,1/5,1/5,0,1/5;
%     1/5,0,1/5,1/5,0,0,1/5,0,1/5,0];
 
 %v = [0;0;1;0;0;1;1;0;0;0];
 %e = [1;1;1;1;1;1;1;1;1;1];
 
 function [y] = PageRank(Q,v)
 
 n = length(v);
 e = ones(n,1);
 
 A = Q + (1/n)*v*transpose(e);
 
 AT = transpose(A);
 
 [V,D] = eig(AT)
 
 end
 
 % find largest magnitude eigenvalue in D, then correspoding eigenvector in
 % V determines PageRank
 