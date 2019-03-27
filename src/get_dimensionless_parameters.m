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

function par=get_dimensionless_parameters(par)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT_PARAMETERS 
%   Input/Output
%       par : array with model parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %----------% Derive other reference parameters %----------%
    par.t0      = par.H*1e3/(par.W0*1e-2);          % time scale [yr]

    %----------% Define other dimensionless parameters %----------%
    par.Deff    = par.D_vol/(1-par.D_vol);              % $D$ : effective partition coeff.
    par.G       = par.Fmax*(par.H/par.Hdry);            % $\Gamma^*$
    par.M       = par.Fmax*(par.H/par.Hdry - 1);        % $\mathcal{M}$
    par.omega   = 2.*pi/par.tp*par.t0;                  % $\omega$
    par.delta0  = par.S0/(par.H*1e3)*par.rhow/par.rhom;
   
    %----------%%----------%%----------%%----------%%----------%
    %----------%  ##  Display input parameters  ##  %----------%
    %----------%%----------%%----------%%----------%%----------%
    if (strcmp(par.verb,'on')==1)
       fprintf('\n\n #-----# \t\t\t    DERIVING DIMENSIONLESS MODEL PARAMETERS   \t\t\t   #-----#') 
       fprintf('\n #-----# Q = %4.2e ; Deff = %6.5e ; Gamma = %4.2f ; M = %4.2f ; omega = %4.2f; delta0 = %4.1e #-----# \n\n',par.Q,par.Deff,par.G,par.M,par.omega,par.delta0); 
    end
    
end