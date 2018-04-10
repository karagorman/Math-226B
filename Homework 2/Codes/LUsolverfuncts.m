% Math 226B - Homework #2
% Problem 4abc

%load('small_ex.mat')
%load('large_ex.mat')

format long e
%load('large_ex1.mat')
%load('large_ex2.mat')
LUfactSolve(2,0,1,2)
 
% [JL,IL,VL] = SparseLowerTri(L);
% c = Lsolve(IL,JL,VL,b);
% [JU,IU,VU] = SparseUpperTri(U);
% x = Usolve(IU, JU, VU, c);


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


% #4
% for part (a)
function [JL,IL,VL] = SparseLowerTri(L)

% using tril
Lt = tril(L,-1);
[JL, KL, VL] = find(Lt);

IL(1) = 1;

for i = 2:size(Lt,2)+1
   countL = nnz(Lt(:,i-1));
   IL(i) = IL(i-1) + countL;
end
IL = transpose(IL);
end


% for part (b)

function [JU,IU,VU] = SparseUpperTri(U)

[JU, KU, VU] = find(U);

IU(1) = 1;

for i = 2:size(U,2)+1
   countU = nnz(U(:,i-1));
   IU(i) = IU(i-1) + countU;
end
IU = transpose(IU);
end

% for part (c)
function c = Lsolve(IL,JL,VL,b)

c = b;

for k = 2:length(IL)-1
    indL = IL(k-1):IL(k) - 1;
    rL = JL(indL);
    c(rL) = c(rL) - VL(indL)*c(k-1);
end

end

% for part (c)
function x = Usolve(IU, JU, VU, c)
x = c;
n = length(c);

for j = length(IU)-1:-1:1
    indU = IU(j):IU(j+1) - 1;
    rU = JU(indU);
    x(j) = x(j)/VU(indU(end));
    x(rU(1:end-1)) = x(rU(1:end-1)) - x(j)*VU(indU(1:end-1));
end

end


