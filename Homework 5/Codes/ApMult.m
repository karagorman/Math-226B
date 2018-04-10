% Math 226B - Homework #5
% Problem 1
% function to compute A'v

function Ap = ApMult(v,gamma,m)

n = m^2;
h = 1/(m+1);

v_bar = A1Mult(v);

v_bar = reshape(v_bar,[m,m]);
Z = fft2DPoisson(m,v_bar);
Z = reshape(Z,[n,1]);

Ap = v + gamma.*Z;


%     function A1v = A1Mult(v)
%         v = reshape(v,[m,m]);
%         Sm = diag(ones(m-1,1),1) - diag(ones(m-1,1),-1);
%         A1v = (h/2)*Sm*v;
%     end

    function z = A1Mult(v)
        z = vertcat(v(m+1:end),zeros(m,1));
        z(m+1:end) = z(m+1:end) - v(1:end-m);
        z = (h/2).*z;
    end


end