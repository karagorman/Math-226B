% Math 226B - Homework #2
% Problem 5
% needs to be run inside the LUsolverFuncts file

function [error,zerosL,zerosU] = LUfactSolve(mat,scale,perm,permType)

if mat == 1
    load('large_ex1.mat')
elseif mat == 2
    load ('large_ex2.mat')
end

if scale == 0 % no scaling
    if perm == 0 % no column reordering
        p0 = 1:size(A,1);
        Ap = A(:,p0);
        
    elseif perm == 1 % yes column reordering
            if permType == 1 % column reordering using colamd
                p0 = colamd(A);
                Ap = A(:,p0);
                
            elseif permType == 2 % column reordering using colperm
                p0 = colperm(A);
                Ap = A(:,p0);  
            end
    end
    [L,U,p,q] = lu(Ap, 'vector');
    n = length(p);
    bp = b(p);

    [JL,IL,VL] = SparseLowerTri(L);
    c = Lsolve(IL,JL,VL,bp);
    [JU,IU,VU] = SparseUpperTri(U);
    d = Usolve(IU, JU, VU, c);

    qi(q) = 1:n;
    p0i(p0) = 1:n;
    x = d(qi);
    x = x(p0i);
    
    % Print results
    zerosL = nnz(L)
    zerosU = nnz(U)
    error = norm(b-A*x)/(norm(b))
    x(10)
    x(100)
    x(1000)
    x(100000)
    x(200000)
    
elseif scale == 1 % yes scaling
    if perm == 0 % no column reordering
        p0 = 1:size(A,1);
        Ap = A(:,p0);
        
    elseif perm == 1 % yes column reordering
            if permType == 1 % column reordering using colamd
                p0 = colamd(A);
                Ap = A(:,p0);
                
            elseif permType == 2 % column reordering using colperm
                p0 = colperm(A);
                Ap = A(:,p0);  
            end
    end
    [L,U,p,q,D] = lu(Ap,'vector');
    n = length(p);
    bp = D\b;
    bp = bp(p);

    [JL,IL,VL] = SparseLowerTri(L);
    c = Lsolve(IL,JL,VL,bp);
    [JU,IU,VU] = SparseUpperTri(U);
    d = Usolve(IU, JU, VU, c);

    qi(q) = 1:n;
    p0i(p0) = 1:n;
    x = d(qi);
    x = x(p0i);

    zerosL = nnz(L)
    zerosU = nnz(U)
    error = norm(b-A*x)/(norm(b))
    x(10)
    x(100)
    x(1000)
    x(100000)
    x(200000)
end
end  

% to compute without using a function.
% #5
% (i) without column reordering, without scaling
% [L,U,p,q] = lu(A,'vector');
% n = length(p);
% bp = b(p);
% 
% [JL,IL,VL] = SparseLowerTri(L);
% c = Lsolve(IL,JL,VL,bp);
% [JU,IU,VU] = SparseUpperTri(U);
% d = Usolve(IU, JU, VU, c);
% 
% qi(q) = 1:n;
% x = d(qi);
% 
% zerosL = nnz(L)
% zerosU = nnz(U)
% error = norm(b-A*x)/(norm(b))
% x(10)
% x(100)
% x(1000)
% x(100000)
% x(200000)


% (i) with scaling
% [L,U,p,q,D] = lu(A, 'vector');
% n = length(p);
% bp = D\b;
% bp = bp(p)
% 
% [JL,IL,VL] = SparseLowerTri(L);
% c = Lsolve(IL,JL,VL,bp);
% [JU,IU,VU] = SparseUpperTri(U);
% d = Usolve(IU, JU, VU, c);
% 
% qi(q) = 1:n;
% x = d(qi);
% 
% zerosL = nnz(L)
% zerosU = nnz(U)
% error = norm(b-A*x)/(norm(b))
% x(10)
% x(100)
% x(1000)
% x(100000)
% x(200000)


% (ii) with column reordering given by p0 = colmad(A), without scaling
% p0 = colamd(A);
% [L,U,p,q] = lu(A(:,p0),'vector');
% bp = b(p);
% 
% 
% [JL,IL,VL] = SparseLowerTri(L);
% c = Lsolve(IL,JL,VL,bp);
% [JU,IU,VU] = SparseUpperTri(U);
% d = Usolve(IU, JU, VU, c);
% 
% qi(q) = 1:n;
% p0i(p0) = 1:n;
% dq = d(qi);
% x = dq(p0i);
% 
% zerosL = nnz(L)
% zerosU = nnz(U)
% error = norm(b-A*x)/(norm(b))
% x(10)
% x(100)
% x(1000)
% x(100000)
% x(200000)

% (ii) with column reordering given by p0 = colmad(A), with scaling
% p0 = colamd(A);
% [L,U,p,q,D] = lu(A(:,p0),'vector');
% n = length(p);
% bp = D\b;
% bp = bp(p);
% 
% [JL,IL,VL] = SparseLowerTri(L);
% c = Lsolve(IL,JL,VL,bp);
% [JU,IU,VU] = SparseUpperTri(U);
% d = Usolve(IU, JU, VU, c);
% 
% qi(q) = 1:n;
% p0i(p0) = 1:n;
% dq = d(qi);
% x = dq(p0i);
% 
% zerosL = nnz(L)
% zerosU = nnz(U)
% error = norm(b-A*x)/(norm(b))
% x(10)
% x(100)
% x(1000)
% x(100000)
% x(200000)


% (iii) with column reordering given by p0 = colperm(A), without scaling
% p0 = colperm(A);
% [L,U,p,q] = lu(A(:,p0),'vector');
% bp = b(p);
% 
% 
% [JL,IL,VL] = SparseLowerTri(L);
% c = Lsolve(IL,JL,VL,bp);
% [JU,IU,VU] = SparseUpperTri(U);
% d = Usolve(IU, JU, VU, c);
% 
% qi(q) = 1:n;
% p0i(p0) = 1:n;
% dq = d(qi);
% x = dq(p0i);
% 
% 3
% zerosL = nnz(L)
% zerosU = nnz(U)
% error = norm(b-A*x)/(norm(b))
% x(10)
% x(100)
% x(1000)
% x(100000)
% x(200000)


% (iii) with column reordering given by p0 = colperm(A), with scaling
% p0 = colperm(A);
% [L,U,p,q,D] = lu(A(:,p0),'vector');
% n = length(p);
% bp = D\b;
% bp = bp(p);
% 
% [JL,IL,VL] = SparseLowerTri(L);
% c = Lsolve(IL,JL,VL,bp);
% [JU,IU,VU] = SparseUpperTri(U);
% d = Usolve(IU, JU, VU, c);
% 
% qi(q) = 1:n;
% p0i(p0) = 1:n;
% dq = d(qi);
% x = dq(p0i);
% 
% zerosL = nnz(L)
% zerosU = nnz(U)
% error = norm(b-A*x)/(norm(b))
% x(10)
% x(100)
% x(1000)
% x(100000)
% x(200000)






