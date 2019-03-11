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

load('column_admittance.mat');

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

mstyle  = {'none','.'};
params_to_plot = [1 , 2];

tpmin = 0.0;
step_marker = 8;

%----------%%----------% Initializing %----------%%----------% 
nfig = nfig + 1; figure(nfig); set(nfig, 'Position', [500, 100, 380, 650]); fig1=nfig;

p = panel(nfig);

% Create panel 2x3
p.pack(2,1);

p.margin=[22 20 5 8]; %% margin[left bottom right top]
p.de.margin = 10;


%----------%%----------% Plotting Admittance  %----------%%----------%
fprintf('\n\n #-----# PLOTTING ADMITANCE AND LAG #-----# \n');

p(1,1).select(); SF_Af=gca;  

for iQ=params_to_plot
    
    fprintf('Plotting admittance for Q=%4.1e \n',par_array(iQ,1).Q);
    
    RelAdm_FLUX = par_array(iQ,1).delta0*abs(FFieldsTop.qh(iQ,:))./MFieldsTop.q(iQ);  
    RelAdm_ECO2 = par_array(iQ,1).delta0*abs(FFieldsTop.qch(iQ,:))./MFieldsTop.qc(iQ);
    
    plot(SF_Af,tperiod(:),RelAdm_FLUX(:)*100,'LineWidth',linew,'color',col_FLUX); hold on;
    plot(SF_Af,tperiod(1:step_marker:end),RelAdm_FLUX(1:step_marker:end)*100,'LineWidth',linew,'LineStyle','none','Marker',mstyle{iQ},'MarkerSize',20,'color',col_FLUX); hold on;
    plot(SF_Af,tperiod(:),RelAdm_ECO2(:)*100,'LineWidth',linew,'color',col_ECO2); hold on;
    plot(SF_Af,tperiod(1:step_marker:end),RelAdm_ECO2(1:step_marker:end)*100,'LineWidth',linew,'LineStyle','none','Marker',mstyle{iQ},'MarkerSize',20,'color',col_ECO2); hold on;
    
    if (iQ==1) 
        hold('on');
        ylabel(SF_Af,{'$A$ [\% per 100-m of SL change]'},'Fontsize',fontsize,'interpreter','latex');
        set(SF_Af,'xlim',[0 200],'ylim',[0 10],'Fontsize',fontsize,'Box','on');
        grid(SF_Af,'on'); SF_Af.XMinorGrid='on'; SF_Af.YMinorGrid='on';
    end
     
end


%----------%%----------% Plotting Lag  %----------%%----------%

p(2,1).select(); SF_Af=gca;

for iQ=params_to_plot

    fprintf('Plotting Lag for Q=%4.1e \n',par_array(iQ,1).Q);

    lag_q  = (wrapTo2Pi(-angle(FFieldsTop.qh(iQ,:))) - pi/2).*tperiod/(2*pi);
    lag_qc = (wrapTo2Pi(-angle(FFieldsTop.qch(iQ,:))) - pi/2).*tperiod/(2*pi);

    if (iQ==1) 
        hold('on');
        xlabel(SF_Af,'$t_p$ [kyr]','Fontsize',fontsize,'interpreter','latex');
        ylabel(SF_Af,{'Lag [kyr]'},'Fontsize',fontsize,'interpreter','latex');
        set(SF_Af,'xlim',[0 200],'ylim',[-5 20],'Fontsize',fontsize,'Box','on');
        grid(SF_Af,'on'); SF_Af.XMinorGrid='on'; SF_Af.YMinorGrid='on';
    end
     
    plot(SF_Af,tperiod(:),lag_q,'LineWidth',linew,'color',col_FLUX); hold on;
    plot(SF_Af,tperiod(1:step_marker:end),lag_q(1:step_marker:end),'LineWidth',linew,'LineStyle','none','Marker',mstyle{iQ},'MarkerSize',20,'color',col_FLUX); hold on;
    plot(SF_Af,tperiod(:),lag_qc,'LineWidth',linew,'color',col_ECO2); hold on;
    plot(SF_Af,tperiod(1:step_marker:end),lag_qc(1:step_marker:end),'LineWidth',linew,'LineStyle','none','Marker',mstyle{iQ},'MarkerSize',20,'color',col_ECO2); hold on;

end
%----------%%----------% Legend %----------%%----------%
h(1)=plot(NaN,NaN,'-','Color',col_FLUX,'linewidth',linew); hold on;
h(2)=plot(NaN,NaN,'-','Color',col_ECO2,'linewidth',linew); hold on;
[hh,icons,plots,txt] = legend(h,{'melt flux','carbon flux'},'Box','off','Fontsize',18,'Position',[0.67 0.85 0.2 0.1],'Units','normalized','Orientation','vertical');  

fprintf('... DONE \n\n')