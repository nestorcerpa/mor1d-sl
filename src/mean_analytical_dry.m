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

function [cs,y,phi,dcs,dy] = mean_analytical_dry(z,par)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MEAN_ANALYTICAL_DRY : Calculatesmean variables and derivatives
%   Inputs
%       z   : depth (can be a scalar or vector)
%       par : array with model parameters
%         par.Deff : effective partition coefficient D=D/(1-D)
%         par.G    : $\Gamma^*$
%         par.M    : $\mathcal{M}$
%         par.Q    : $\mathcal{Q}$
%    Outputs
%       cs  : mean solid concentration
%       y   : change of variable y = Q(phi)
%       phi : mean porosity
%       dcs  : spatial-derivative of cs
%       dy  : spatial-derivative of y
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    D = par.Deff;
    G = par.G;
    M = par.M;
    n = par.n;
    Q = par.Q;

    if M>1e-6
        error('Error: not dry model, check par.M')
    end

    %----------% Calculate base state functions
    y=G*z;

    %----------% Determine phi (requires solution of an algebraic equation)
    phi=zeros(size(y));
    for yt=y
        if n==2 %special case can be done with roots function
            pol_p = [Q, -2*Q, Q, 1,-yt];
            r_p=roots(pol_p);
            r_p = real(r_p(abs(imag(r_p))<1e-14));
            r_p=r_p(r_p<=1 & r_p>=0);
            if numel(r_p)~=1
                disp(r_p)
                disp(pol_p)
                pause(0.1)
                error('incorrect number of roots found in the range (0,1)')
            end
            phi(yt==y)=r_p;
        elseif n==3 %special case can be done with roots function
            pol_p = [Q, -2*Q, Q, 0, 1,-yt];
            r_p=roots(pol_p);
            r_p = real(r_p(abs(imag(r_p))<1e-14));
            r_p=r_p(r_p<=1 & r_p>=0);
            if numel(r_p)~=1
                disp(r_p)
                error('incorrect number of roots found in the range (0,1)')
            end
            phi(yt==y)=r_p;
        else %fallback for n not integer
            phi(yt==y)=fzero(@(p) Q*p^n*(1-p)^2+p -yt,[0 1]);
        end
    end

    %----------% Calculate concentration (change variables)
    cs=D./(D+y);

    %----------% Calculate derivatives
    dcs=-cs*G./(D+G*z);
    dy=G;


end
