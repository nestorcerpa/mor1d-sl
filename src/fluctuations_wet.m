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

function [xh,yh,csh,phih] = fluctuations_wet(z_out,par)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FLUCTUATIONS Calculates flucuating variables 
%   Inputs
%       z_out is set of depths at which outputs are evaluated (can be a scalar or vector)
%       par : array with model parameters
%         par.Deff : effective partition coefficient D=D/(1-D)
%         par.G    : $\Gamma^*$
%         par.M    : $\mathcal{M}$
%         par.Q    : $\mathcal{Q}$
%   Outputs
%       xh   : $\hat{x}$=M/c;
%       yh   : fluctuating liquid flux $\hat{y}=\hat{Q}$
%       csh  : fluctuating carbon concentration in solid \hat{c_s}
%       phih : fluctuating porosity \hat{\phi} 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if strcmp(par.verb,'on')==1
       fprintf('\n Solving ODEs for fluctuations ...\n') 
    end
    
    D = par.Deff; 
    G = par.G;
    M = par.M;
    n = par.n;
    Q = par.Q;
    w = par.omega;

    %----------% specifiy ode45 options (relatively tight tolerances)
    options = odeset('AbsTol',par.ODEopts.AbsTol,'RelTol',par.ODEopts.RelTol,'Stats',par.verb);

    %----------% calculate mean state derivatives (needed for boundary conditions)
    z=0;
    [~,~,~,~,dx0,dy0] = mean_analytical_wet(z,par);
    
    %----------% calculate fluctuating state
    sol = ode45(@(zi,y) ode(zi,y,par),[0 1],[-dx0;-dy0],options);

    %----------% output
    xh=deval(z_out,sol,1);
    yh=deval(z_out,sol,2);
    [~,~,~,phi,~,~] = mean_analytical_wet(z_out,par);
    dQdphi =  Q*(n*phi.^(n-1).*(1-phi).^2 - 2*(1-phi).*phi.^n)+1; 
    phih=yh./dQdphi;
    csh = xh/M; 
    
    if (strcmp(par.verb,'on')==1); fprintf('... DONE \n\n'); end;
    
        function deriv = ode(z,y,par)
            xh = y(1);
            yh = y(2);
            [x,y,~,phi,dx,dy] = mean_analytical_wet(z,par);
            dQdphi =  Q*(n*phi^(n-1)*(1-phi)^2 - 2*(1-phi)*phi^n)+1;
            phih=yh/dQdphi;
             if (strcmp(par.Gammap,'on')==1)
                dxh = (D+x+y)^(-1) *( 1i*w*G*x - xh*(1i*w*(D+phi+x)+dy) - yh*dx);
                dyh = 1i*w*(xh-G-phih) +dxh;
            elseif (strcmp(par.Gammap,'off')==1)
                dxh = (D+y)^(-1) *(- xh*(1i*w*(D+phi+x)+dy) - yh*dx);
                dyh = 1i*w*(xh-phih) +dxh;         
            end
            deriv = [dxh;dyh];
        end

 
end

