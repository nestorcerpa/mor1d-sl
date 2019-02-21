function Y = idft(F)
% IDFT   INVERSE DISCRETE FOURIER TRANSFORM
% 
% IDFT(F) takes in the ouput of the discrete Fourier transform function dft
% and reconstructs a time series Y. Input F is a structure with the form:
%   F.gamma  = coefficients of complex exponential terms for k=1:(N-1)/2
%   F.alpha0 = mean of the time-series (alpha_0)

N = length(F.gamma)*2+1;
n = 0:N-1;
Y = F.alpha0*ones(1,N);

for k = 1:length(F.gamma)
    Y = Y + F.gamma(k)*exp(1i*2*pi*k*n/N);
end
Y = real(Y);