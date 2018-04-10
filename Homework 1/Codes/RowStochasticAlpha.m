function [yalpha] = RowStochasticAlpha(n,Q,Jv,alpha,x0)

e = ones(n,1);
v = zeros(n,1);

for i=1:length(Jv)
    v(Jv(i)) = 1;
end

v = sparse(v);
vT = transpose(v);
vTx0 = vT*x0;
QT = transpose(Q);
eTx0 = transpose(e)*x0;

ATx0 = QT*x0 + (1/n)*e*vTx0;
ATx0_alpha = alpha*ATx0 + (1 - alpha)*(1/n)*e*eTx0;

yalpha = ATx0_alpha;
end

