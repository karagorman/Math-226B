function [y_fast] = fftCirculant(c,x)

n = length(x);

%C = [2,-1,0,-3;
 %    -3,2,-1,0;
  %   0,-3,2,-1;
  %   -1,0,-3,2];
 
 %c = C(1,:);
 
 lambdaVec = conj(fft(c));
 
 lambdaMat = zeros(n,n);
 
 for i=1:n
     lambdaMat(i,i) = lambdaVec(i);
 end
 
 
 %x = [-1;2;1;4];
 
 %y_exact = C*x;
 y_fast = (1/n)*conj(fft(conj(lambdaVec.*fft(x))));
end
 