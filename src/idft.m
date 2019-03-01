% % % MIT Licence
% % % 
% % % Copyright (c) 2019
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

end