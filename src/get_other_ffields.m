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

function [Wh,fh,fch] = get_other_ffields(phib,csb,phih,csh,par)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GET_OTHER_FFIELDS Calculates base state variables and derivative
%   Inputs
%       par : array with model parameters
%             par.Q : \mathcal{Q} 
%             par.n : n
%       csb  : mean carbon concentration in solid \bar{c}_s
%       phib : mean porosity \bar{\phi}
%       csh  : fluctuating carbon concentration in solid \hat{c_s}
%       phih : fluctuating porosity \hat{\phi} 
%   Outputs
%       Wh   : fluctuating solid velocity \hat{W}
%       fh   : fluctuating melt flux \hat{f}
%       fch  : fluctuating carbon flux \hat{f}_c
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Q = par.Q;
    n = par.n;

    %----------% Fluctuating Solid velocity %----------%
    %   \hat{W} = \cal{Q} (\bar{[phi}^n \hat{\phi} - n \hat{\phi}\bar{\phi}^(n-1)*(1-\bar{\phi}))
    %
    Wh = Q * (phib.^n .* phih - n * phih .* phib.^(n-1) .* (1-phib));
    
    %----------% Fluctuating melt flux %----------% 
    Wb = 1 - Q .* phib.^n.*(1-phib); % \bar{W}
    fh = -(1-phib).*Wh + Wb .* phih;

    %----------% Fluctuating chemical flux %----------% 
    fb = 1 - (1-phib).*Wb;
    fch = (csb.*fh + fb.*csh)/par.D_vol;
    
end