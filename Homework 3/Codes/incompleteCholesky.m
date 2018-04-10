% Math 226B - HW #3
% Problem 5b
% Write a Matlab function that efficiently computes the entries 
% ljk, (j,k)?E, (8)
% of the incomplete Cholesky factor L of A. The input of your function 
% should be A and the integer vectors J and I from (a), and the output 
% should be a vector VL that contains the entries (8) in the column-wise 
% order given by J and I.


function [VL] = incompleteCholesky(J,I,VA,A)

VL = VA;

for k=1:size(A,1)
    VL(I(k)) = sqrt(VL(I(k)));
    ind = I(k):I(k+1)-1;
    ind = ind(2:end);
    rowInd = J(ind);
    
    VL(ind) = VL(ind)/VL(I(k));
    
    for j = 1:length(ind)
        i = rowInd(j);
        indInd = I(i):I(i+1) - 1;
        rowIndInd = J(indInd);
        [~,j1,j2] = intersect(rowInd, rowIndInd); % idea from Karry
        VL(indInd(j2)) = VL(indInd(j2)) - VL(ind(j1)).*VL(ind(j));
    end
end
end
    
     
    
    
    
