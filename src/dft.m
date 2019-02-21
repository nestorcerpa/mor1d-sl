function F = dft(y)
% DFT   Discrete Fourier transform
%   DFT(Y) computes the Discrete Fourier transform of an input
%   time-series Y.  F=DFT(Y) uses the built-in FFT and returns a
%   structure F as follows:
%     F.alpha0 = mean of the time-series
%     F.alpha  = coefficients of cosine terms for k=1:Nyquist
%     F.beta   = coefficients of the sine terms for k=1:Nyquist
%     F.power  = normalised power-spectrum of the time-series

  %1. Compute the discrete Fourier transform using Matlab's fft
  G = fft(y);

  %2. Determine and select the unique entries
  N  = length(y);
  Nu = ceil((N+1)/2);
  G  = G(1:Nu);
  
  %3. Scale value appropriately and assemble output structure
  F.alpha0 = real(G(1))/N;
  F.alpha  =  2*real(G(2:end))/N;
  F.beta   = -2*imag(G(2:end))/N;
  F.gamma  =  2*G(2:end)/N;
  F.power  = 0.5*(F.alpha.^2 + F.beta.^2)/var(y); 
