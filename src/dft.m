% % % MIT Licence
% % % 
% % % Copyright (c) 2019
% % %     Nestor G. Cerpa       (University of Montpellier) [nestor.cerpa@gm.univ-montp2.fr]
% % %     David W. Rees Jones   (University of Oxford)      [david.reesjones@earth.ox.ac.uk]
% % %     Richard F. Katz       (University of Oxford)      [richard.katz@earth.ox.ac.uk] 
% % % 
% % % Permission is hereby granted, free of charge, to any person obtaining a copy
% % % of this software and associated documentation files (the "Software"), to deal
% % % in the Software without restriction, including without limitation the rights
% % % to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% % % copies of the Software, and to permit persons to whom the Software is
% % % furnished to do so, subject to the following conditions:
% % % 
% % % The above copyright notice and this permission notice shall be included in all
% % % copies or substantial portions of the Software.
% % % 
% % % THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% % % IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% % % FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% % % AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% % % LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% % % OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% % % SOFTWARE.

function F = dft(y)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DFT   Discrete Fourier transform
%   DFT(Y) computes the Discrete Fourier transform of an input
%   time-series Y.  F=DFT(Y) uses the built-in FFT and returns a
%   structure F as follows:
%     F.alpha0 = mean of the time-series
%     F.alpha  = coefficients of cosine terms for k=1:Nyquist
%     F.beta   = coefficients of the sine terms for k=1:Nyquist
%     F.power  = normalised power-spectrum of the time-series
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
end