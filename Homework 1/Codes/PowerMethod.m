% Math 226B - Homework #1
% Problem 4b
% Write a Matlab program that implements the power method, as stated in the
% notes provided on the course website.


function [y] = PowerMethod(x, epsilon, kMax, n, Q, Jv, alpha)

lambda = [];
x_abs = abs(x);
[~,b] = max(x_abs);
X = x(b);
lambda(1) = X;

itCount = 0;

%%%
e = ones(n,1);
v = zeros(n,1);

for i=1:length(Jv)
    v(Jv(i)) = 1;
end

v = sparse(v);
vT = transpose(v);
QT = transpose(Q);
%%%

for k=2:kMax
    vTx = vT*x;
    ATx = QT*x + (1/n)*e*vTx;
    eTx = transpose(e)*x;
    x = alpha*ATx + (1-alpha)*(1/n)*e*eTx;
    x_abs = abs(x);
    %[~,b] = sort(x_abs,'descend');
    [~,b] = max(x_abs);
    X = x(b);
    lambda(k) = X;
    x = x/lambda(k);
    itCount = itCount + 1;
    
    if (abs(lambda(k) - lambda(k-1)) <= epsilon*(abs(lambda(k-1))))
        break
    end
end

itCount = itCount
evalApprox = lambda(itCount+1)
evecApprox = x;
vTx = vT*x;
ATx = QT*x + (1/n)*e*vTx;
eTx = transpose(e)*x;
ATalphax = alpha*ATx + (1-alpha)*(1/n)*e*eTx;
relEigResidual = norm((ATalphax - evalApprox*x),inf)/norm(x,inf)
v = x;
[~,d] = sort(v, 'descend');
pageRank = d(1:10)

end
