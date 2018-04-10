% Math 226B - Homework #2
% Problem 4a

%load('small_ex.mat')
load('large_ex.mat')

%[J,I,VL] = SparseLowerTri(L);
[J,I,VU] = SparseUpperTri(U);


% for part (a)
function [J,I,VL] = SparseLowerTri(L)

[row, col, v] = find(L);

J = row;
VL = v;
K = col;

I = find(J-K == 0);
i = length(I);
I(i+1) = length(J) + 1;
end


% for part (b)

function [J,I,VU] = SparseUpperTri(U)

[row, col, v] = find(U);

J = row;
VU = v;
I(1) = 1;

for i = 2:size(U,2)+1
   count = nnz(U(:,i-1));
   I(i) = I(i-1) + count;
end
end




