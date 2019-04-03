% MIT Licence
% 
% Copyright (c) 2019
%     Nestor G. Cerpa       (University of Montpellier) [nestor.cerpa@gm.univ-montp2.fr]
%     David W. Rees Jones   (University of Oxford)      [david.reesjones@earth.ox.ac.uk]
%     Richard F. Katz       (University of Oxford)      [richard.katz@earth.ox.ac.uk] 
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.

close all; clear all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----------% Reading data from  file %----------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%----------% Loading files %----------%
load('mor1d_admittance.mat');

par_array  = data.par_array;
MFieldsTop = data.MFieldsTop; 
FFieldsTop = data.FFieldsTop; 
tperiod    = data.tp_array;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----------%%----------% FIGURE %----------%%----------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nfig = 0; 

%----------%%----------% Options %----------%%----------%

linew = 3;      % Line width for plots
legendsize= 16; % Legend fontsize
fontsize = 16;  % Labels fontsize
marksize = 20;

col_FLUX = [0 0 1];
col_ECO2 = [0.5 0 0.5];

params_to_plot = [1 ];

tpmin = 0.0;
step_marker = 8;

%----------%%----------% Initializing %----------%%----------% 
nfig = nfig + 1; figure(nfig); set(nfig, 'Position', [200, 100, 650, 350]); 

p = panel(nfig);

% Create panel 2x3
p.pack(1,2);

p.margin=[18 20 3 8]; %% margin[left bottom right top]
p.de.margin = 18;


%----------%%----------% Plotting Admittance  %----------%%----------%
fprintf('\n\n #-----# PLOTTING ADMITANCE AND LAG #-----# \n');

p(1,1).select(); SF=gca;  

for iQ=params_to_plot
    
    fprintf('Plotting admittance for Q=%4.1e \n',par_array(iQ,1).Q);
    
    %----------% Parameters %----------%
    par=data.par_array(iQ);
    n = par.n;
    D = par.Deff;
    phibar = MFieldsTop.phi(iQ,end);
    cbar   = MFieldsTop.cs(iQ,end)*par.M;
    Qbar   = MFieldsTop.q(iQ,end);     
    dQ_dphi = par.Q*(n*phibar^(n-1)*(1-phibar)^2 - 2*phibar^n*(1-phibar) ) + 1.0;    
    w0 = par.W0*par.Fmax/((par.Fmax/par.Q)^(1/par.n));  % in cm/yr

    %----------% Admittance analytic %----------%
    % Admittance approximation
    Afactor = (D+phibar)/(D+phibar+cbar);
    Af0 = par.delta0*par.G/Qbar * (Afactor) * dQ_dphi;                          % Eq. (C.5)        
    Afc0 = par.delta0*par.G*(Afactor*dQ_dphi./Qbar - 1./(D+phibar+cbar) );      % Eq. (C.7)
    % Admittance approximation for D<<\phi<<1
    Af_smallD  = n*par.delta0*par.Fmax^(1-1/n)*par.H/par.Hdry*par.Q^(1/n);      % Eq. (C.6)
    Afc_smallD = (n-1)*par.delta0*par.Fmax^(1-1/n)*par.H/par.Hdry*par.Q^(1/n);  % Eq. (C.8)
    % Phase approximation
    Imf  = (n-1)*par.G*(par.Q/par.Fmax)^(1/n);
    Imfc = ((n^2-2*n+1/n)/(n-1))*par.G*(par.Q/par.Fmax)^(1/n);
    omega= (2*pi*par.t0)./(tperiod*1e3);

    % Admittance from models
    RelAdm_FLUX = par_array(iQ,1).delta0*abs(FFieldsTop.qh(iQ,:))./MFieldsTop.q(iQ);  
    RelAdm_ECO2 = par_array(iQ,1).delta0*abs(FFieldsTop.qch(iQ,:))./MFieldsTop.qc(iQ);
 
    %----------% Critical period %----------%
    tp_star = 1/(n-1)*par.Hdry/(w0*1e-5)*1e-3;      % in kyr
    plot(SF,[tp_star,tp_star],get(SF,'ylim'),'LineWidth',linew-1,'color','k','LineStyle','--'); hold on;
    
    %----------% Plotting admittance %----------%
    % Melt flux
    plot(SF,tperiod,2*RelAdm_FLUX(:)*100,'LineWidth',linew,'color',col_FLUX); hold on;
    plot(SF,tperiod,2*Af0*100+zeros(size(tperiod)),'LineWidth',linew-1,'color',col_FLUX,'LineStyle',':'); hold on;          % (C.5) 
    plot(SF,tperiod,2*Afc0*100+zeros(size(tperiod)),'LineWidth',linew-1,'color',col_ECO2,'LineStyle',':'); hold on;         % (C.7)
    % Carbon flux
    plot(SF,tperiod,2*RelAdm_ECO2(:)*100,'LineWidth',linew,'color',col_ECO2); hold on;
    plot(SF,tperiod,2*Af_smallD*100+zeros(size(tperiod)),'LineWidth',linew-1,'color',col_FLUX,'LineStyle','--'); hold on;   % (C.6)
    plot(SF,tperiod,2*Afc_smallD*100+zeros(size(tperiod)),'LineWidth',linew-1,'color',col_ECO2,'LineStyle','--'); hold on;  % (C.8)
    
    if (iQ==1) 
        
        %----------% Axis appearance %----------%
        hold('on');
        xlabel(SF,'$t_p$ [kyr]','Fontsize',fontsize,'interpreter','latex');
        ylabel(SF,{'$A$ [\% per 100-m of SL change]'},'Fontsize',fontsize,'interpreter','latex');
        set(SF,'XScale','log','YScale','log','xlim',[1 200],'ylim',[1 30],'Fontsize',fontsize,'Box','on');
        grid(SF,'on'); SF.XMinorGrid='on'; SF.YMinorGrid='on';
        
        %----------% Critical period %----------%
        tp_star = 1/(par.n-1)*par.Hdry/(w0*1e-5)*1e-3;      % in kyr
        plot(SF,[tp_star,tp_star],get(SF,'ylim'),'LineWidth',linew-1,'color','k','LineStyle','--'); hold on;
        
    end
     
end


%----------%%----------% Plotting Lag  %----------%%----------%

p(1,2).select(); SF=gca;

for iQ=params_to_plot
    
    fprintf('Plotting Lag for Q=%4.1e \n',par_array(iQ,1).Q);

    %----------% Parameters %----------%
    w0 = par.W0*par.Fmax/((par.Fmax/par.Q)^(1/par.n));  % in cm/yr

    
    lag_q_asym  = (wrapTo2Pi(-atan(Imf./omega)-pi) - pi/2).*tperiod/(2*pi);
    lag_qc_asym  = (wrapTo2Pi(-atan(Imfc./omega)-pi) - pi/2).*tperiod/(2*pi);
    lag_q  = (wrapTo2Pi(-angle(FFieldsTop.qh(iQ,:))) - pi/2).*tperiod/(2*pi);
    lag_qc = (wrapTo2Pi(-angle(FFieldsTop.qch(iQ,:))) - pi/2).*tperiod/(2*pi);
        
    %----------% Plotting lag %----------%
    plot(SF,tperiod,lag_q,'LineWidth',linew,'color',col_FLUX); hold on;
    plot(SF,tperiod,lag_qc,'LineWidth',linew,'color',col_ECO2); hold on;
    plot(SF,tperiod,lag_q_asym,'LineWidth',linew-1,'color',col_FLUX,'LineStyle','--'); hold on;
    plot(SF,tperiod,lag_qc_asym,'LineWidth',linew-1,'color',col_ECO2,'LineStyle','--'); hold on;

    if (iQ==1) 
        
        %----------% Axis appearance %----------%
        hold('on');
        xlabel(SF,'$t_p$ [kyr]','Fontsize',fontsize,'interpreter','latex');
        ylabel(SF,{'Lag [kyr]'},'Fontsize',fontsize,'interpreter','latex');
        set(SF,'XScale','log','xlim',[1 200],'ylim',[-5 20],'Fontsize',fontsize,'Box','on');
        grid(SF,'on'); SF.XMinorGrid='on'; SF.YMinorGrid='on';
        
        %----------% Critical period %----------%
        tp_star = 1/(par.n-1)*par.Hdry/(w0*1e-5)*1e-3;      % in kyr
        plot(SF,[tp_star,tp_star],get(SF,'ylim'),'LineWidth',linew-1,'color','k','LineStyle','--'); hold on;
        
    end
    
end
%----------%%----------% Legend %----------%%----------%
h(1)=plot(NaN,NaN,'-','Color',col_FLUX,'linewidth',linew); hold on;
h(2)=plot(NaN,NaN,'-','Color',col_ECO2,'linewidth',linew); hold on;
[hh,icons,plots,txt] = legend(h,{'melt flux','carbon flux'},'Box','off','Fontsize',18,'Location','southwest','Units','normalized','Orientation','vertical');  

fprintf('... DONE \n\n')
