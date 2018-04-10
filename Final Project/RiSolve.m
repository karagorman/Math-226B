
% solve on Ri
% inputs: kmin, kmax, s, c,ic,d,id

function w = RiSolve(Ri_xind,Ri_yind,kmin,kmax,s,c,ic,d,id,Num_Mat,m)

W = zeros(size(Num_Mat));
% reshape to a matrix
for k = 1:(kmax - kmin + 1)
    [data_x,data_y] = find(Num_Mat == (kmin + k - 1));
    W(data_x,data_y) = s(k);
end

% put into fft solver
W(Ri_xind,Ri_yind) = fft2DPoisson(m,W(Ri_xind,Ri_yind),c,ic,d,id);

w = zeros(length(s),1);
% reshape into vector
for k = 1:(kmax - kmin + 1)
    [data_x,data_y] = find(Num_Mat == (kmin + k - 1));
    w(k) = W(data_x,data_y);
end

end
    
    