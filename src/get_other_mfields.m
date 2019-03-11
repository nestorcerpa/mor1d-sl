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

function [W,q,qc] = get_other_mfields(phi,cs,par)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GET_OTHER_MFIELDS Calculates base state variables and derivative
%   Inputs
%       par : array with model parameters
%             par.Q : \mathcal{Q} 
%             par.n : n
%       cs  : mean carbon concentration in solid \bar{c}_s
%       phi : mean porosity \bar{\phi}
%   Outputs
%       W   : mean solid velocity \bar{W}
%       q   : mean melt flux \bar{q}
%       qc  : mean carbon flux \bar{q}_c
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Q = par.Q;
    n = par.n;
    
    %----------% Mean solid velocity %----------% 
    W  = 1 - Q .* (phi.^n.*(1-phi)); % \bar{W}
    
    %----------% Mean melt flux %----------% 
    q  = 1 - (1-phi).*W;
    
    %----------% Mean chemical flux %----------% 
    qc = q.*cs/par.D_vol;
    
end
