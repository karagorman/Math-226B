% Math 226B - Homework #1
% Problem 4a
% Write two Matlab functions, which use n, Q, and Jv as inputs, to compute 
% matrix- vector products of the form: 
% 1. y = AT x, where x ? Rn is a dense vector in general, and
% 2. y = AT? x, where x ? Rn is a dense vector in general,
% respectively, as efficiently as possible. Here, A?, 0 ? ? ? 1, is the 
% family of matrices defined in (1)
%
% Use your Matlab functions to compute the matrix-vector products:
% y = AT e, y0.5 = AT0.5e, and y0.85 = AT0.85e
% for the case of n = 10 websites stated in Problem 2, and print out y, 
% y0.5, and y0.85.
% 
% Next, run your Matlab functions on the large graph G with n = 685230 
% nodes and 7600595 edges given in ?www0.mat?. For the vector x0 given in 
% ?x0.mat?, compute the matrix-vector products
% y = AT x0 and y0.85 = AT0.85x0,
% and print out the entries with indices 2, 222222, 300000, and 400000 of 
% both vectors y and y0.85.

function [y] = RowStochastic(n,Q,Jv,x0)

e = ones(n,1);
v = zeros(n,1);

for i=1:length(Jv)
    v(Jv(i)) = 1;
end

v = sparse(v);
vT = transpose(v);
vTx0 = vT*x0;
QT = transpose(Q);

ATx0 = QT*x0 + (1/n)*e*vTx0;

y = ATx0;
end
