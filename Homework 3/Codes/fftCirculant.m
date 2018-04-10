function [y_fast] = fftCirculant(c,x)

 n = length(x);
 
 lambdaVec = conj(fft(c'));
 
 y_fast = conj(fft(conj(lambdaVec.*fft(x))))/n;
 
end
 