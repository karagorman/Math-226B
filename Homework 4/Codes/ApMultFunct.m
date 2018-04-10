function Avp = ApMultFunct(L,U,D,D1,z)

us = U\z;
d = diag(D1).*us;
v = z + d;
ls = L\v;
sum = ls + us;
Avp = D*sum;
end