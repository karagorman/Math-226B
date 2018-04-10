% Math 226B - Homework #2
% Problem 1
% For m = m0 := 67, compute the Cholesky factor L of A for the following 
% 5 cases:
% For all 5 runs, report the number of nonzero entries of the Cholesky 
% factor L and submit a plot of the sparsity structure of L (as produced 
% by the Matlab command spy(L)).

%A = make_2d_laplacian(67);
A = make_3d_laplacian(37);

hold on
figure
% (i) No reordering of the rows and columns of A
L1 = chol(A,'lower');
non_zeros1 = nnz(L1)
subplot(2,3,1)
spy(L1)
title('No Reordering')

% (ii) Reordering with symamd;
p2 = symamd(A);
L2 = chol(A(p2,p2),'lower');
non_zeros2 = nnz(L2)
subplot(2,3,2)
spy(L2)
title('Reordering with symamd')


% (iii) Reordering with colamd; 
p3 = colamd(A);
L3 = chol(A(p3,p3),'lower');
non_zeros3 = nnz(L3)
subplot(2,3,3)
spy(L3)
title('Reordering with colamd')


% (iv) Reordering with symrcm;
p4 = symrcm(A);
L4 = chol(A(p4,p4),'lower');
non_zeros4 = nnz(L4)
subplot(2,3,4)
spy(L4)
title('Reordering with symrcm')


% (v) Reordering with colperm.
p5 = colperm(A);
L5 = chol(A(p5,p5),'lower');
non_zeros5 = nnz(L5)
subplot(2,3,5)
spy(L5)
title('Reordering with colperm')


% For m = (2^i)*m0, i = 0,1,..., compute the Cholesky factor L of A using the 
% symamd ordering. As i increases, you will run out of storage or CPU time. 
% For which i does this happen on the machine you run Matlab on?


%m0 = 67;
m0 = 37;

for i = 0:10
    i
    m = (2^i)*m0;
    %A = make_2d_laplacian(m);
    A = make_3d_laplacian(m);
    p = symamd(A);
    L = chol(A(p,p),'lower');
end

%m0 = 67;
m0 = 37;
i = 0;

while (i >= 0)
    m = (2^i)*m0;
    i = i + 1
    %A = make_2d_laplacian(m);
    A = make_3d_laplacian(m);
    p = symamd(A);
    L = chol(A(p,p),'lower');
end

