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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----------%%----------% Reading data from  file %----------%%----------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%----------% Loading files %----------%
load('mor1d_admittance.mat');

par_array  = data.par_array;
BotSurfLag = data.lagBCtoSurf;
tperiod    = data.tp_array;


prompt = '\n\n Admittance to plot \n 1: Wet \n 2: Dry \n 3: Dry basal flux \n? ';
opt_admit = input(prompt);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----------%%----------% FIGURE %----------%%----------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nfig = 0; 

%----------%%----------% Options %----------%%----------%

linew = 3;      % Line width for plots
legendsize= 16; % Legend fontsize
fontsize = 16;  % Labels fontsize
marksize = 20;

col_FLUX   = [0.0 0.0 1.0];
col_ECO2   = [0.5 0.0 0.5];
col_RMOR0  = [0.4 0.6 1.0];
col_RMOR1  = [0.0 0.2 0.4];
col_RCMOR0 = [1.0 0.6 1.0];
col_RCMOR1 = [0.8 0.0 0.4];

switch opt_admit
    case 1 
        params_to_plot = [1 , 3];
    case 2 
        params_to_plot = [4];
    case 3 
        params_to_plot = [5];
end
lins    = {'-';   ':';   '--'; '-'; '-'};

tpmin = 0.0;
step_marker = 8;

%----------%%----------% Initializing %----------%%----------% 
nfig = nfig + 1; figure(nfig); set(nfig, 'Position', [500, 100, 650, 350]); fig1=nfig;

p = panel(nfig);

% Create panel 2x3
p.pack(1,2);

p.margin=[22 20 5 6]; %% margin[left bottom right top]
p.de.margin = 18;


%----------%%----------% Plotting Admittance  %----------%%----------%
fprintf('\n\n #-----# PLOTTING ADMITANCE AND LAG #-----# \n');

p(1,1).select(); SF=gca;  

for iQ=params_to_plot
    
    fprintf('Plotting admittance for Q=%4.1e \n',data.par_array(iQ,1).Q);
    
    %%% Fluctuations to plot
    RelAdm_FLUX      = data.par_array(iQ,1).delta0*abs(data.FFieldsTop.qh(iQ,:))./data.MFieldsTop.q(iQ);  
    RelAdm_ECO2      = data.par_array(iQ,1).delta0*abs(data.FFieldsTop.qch(iQ,:))./data.MFieldsTop.qc(iQ);
    RelAdm_RMOR0  = data.par_array(iQ,1).delta0*abs(data.FFieldsMOR.Rmorh0(iQ,:))./data.MFieldsMOR.Rmor(iQ);
    RelAdm_RMOR1  = data.par_array(iQ,1).delta0*abs(data.FFieldsMOR.Rmorh(iQ,:))./data.MFieldsMOR.Rmor(iQ);  
    RelAdm_RCMOR0 = data.par_array(iQ,1).delta0*abs(data.FFieldsMOR.Rcmorh0(iQ,:))./data.MFieldsMOR.Rcmor(iQ);
    RelAdm_RCMOR1 = data.par_array(iQ,1).delta0*abs(data.FFieldsMOR.Rcmorh(iQ,:))./data.MFieldsMOR.Rcmor(iQ);     
    
    if (opt_admit ~= 3)  % If wet or dry models
        if (opt_admit == 1) 
            plot(SF,tperiod,2*RelAdm_FLUX*100,'LineWidth',linew,'color',col_FLUX,'LineStyle',lins{iQ}); hold on;
        elseif (opt_admit == 2)
            RelAdm_FLUX_BK15 = data.par_array_BK15(1,1).delta0*abs(data.FFieldsTop_BK15.qh(1,:))./data.MFieldsTop_BK15.q(1);
            plot(SF,tperiod,2*RelAdm_RMOR0*100,'LineWidth',linew,'color',col_RMOR0,'LineStyle',lins{iQ}); hold on;
            plot(SF,tperiod,2*RelAdm_RMOR1*100,'LineWidth',linew,'color',col_RMOR1,'LineStyle',lins{iQ}); hold on;
            plot(SF,tperiod,2*RelAdm_FLUX*100,'LineWidth',linew,'color',col_FLUX,'LineStyle',lins{iQ}); hold on;
        end
    end
    if  (opt_admit ~= 2) % If wet or basal-flux models
        tpmin=8.0; [~,idx_min]  = min(abs(tperiod(:)-tpmin));
        if (opt_admit == 1) 
            plot(SF,tperiod(idx_min:end),2*RelAdm_ECO2(idx_min:end)*100,'LineWidth',linew,'color',col_ECO2,'LineStyle',lins{iQ}); hold on;
        elseif (opt_admit == 3) 
            RelAdm_ECO2_BK15 = data.par_array_BK15(1,1).delta0*abs(data.FFieldsTop_BK15.qch(2,:))./data.MFieldsTop_BK15.qc(2);
            plot(SF,tperiod(idx_min:end),2*RelAdm_ECO2(idx_min:end)*100,'LineWidth',linew,'color',col_ECO2,'LineStyle',lins{iQ}); hold on;    
            plot(SF,tperiod(idx_min:end),2*RelAdm_ECO2_BK15(idx_min:end)*100,'LineWidth',linew-2,'color',col_ECO2,'LineStyle',lins{iQ}); hold on; 
        end
    end
    
    if (iQ==params_to_plot(1)) 
        hold('on');
        xlabel(SF,'$t_p$ [kyr]','Fontsize',fontsize,'interpreter','latex')
        ylabel(SF,{'$A$ [\% per 100-m of SL change]'},'Fontsize',fontsize,'interpreter','latex');
        set(SF,'xlim',[0 200],'Fontsize',fontsize,'Box','on');
         if (opt_admit ~= 3)
            set(SF,'ylim',[0 20]);
        else
            set(SF,'ylim',[0 80]);
        end
        grid(SF,'on'); SF.XMinorGrid='on'; SF.YMinorGrid='on';
    end
     
end


%----------%%----------% Plotting Lag  %----------%%----------%

p(1,2).select(); SF=gca;

for iQ=params_to_plot

    fprintf('Plotting Lag for Q=%4.1e \n',par_array(iQ,1).Q);

    lag_q     = (wrapTo2Pi(-angle(data.FFieldsTop.qh(iQ,:)))  - pi/2).*tperiod/(2*pi);
    lag_qc    = (wrapTo2Pi(-angle(data.FFieldsTop.qch(iQ,:))) - pi/2).*tperiod/(2*pi);
    lag_Rmor0 = (wrapTo2Pi(-angle(data.FFieldsMOR.Rmorh0(iQ,:))) - pi/2).*tperiod/(2*pi);
    lag_Rmor1 = (wrapTo2Pi(-angle(data.FFieldsMOR.Rmorh(iQ,:)))  - pi/2).*tperiod/(2*pi);
    
    if (opt_admit ~= 3)  % If wet or dry models
        if (opt_admit == 1) 
            plot(SF,tperiod(:),lag_q,'LineWidth',linew,'color',col_FLUX,'LineStyle',lins{iQ}); hold on;
        elseif (opt_admit == 2)
            lag_q_BK15   = (wrapTo2Pi(-angle(data.FFieldsTop_BK15.qh(1,:))) - pi/2).*tperiod/(2*pi);
            plot(SF,tperiod,lag_Rmor0,'LineWidth',linew,'color',col_RMOR0,'LineStyle',lins{iQ}); hold on;    
            plot(SF,tperiod,lag_Rmor1,'LineWidth',linew,'color',col_RMOR1,'LineStyle',lins{iQ}); hold on;   
            plot(SF,tperiod,lag_q,'LineWidth',linew,'color',col_FLUX,'LineStyle',lins{iQ}); hold on; 
            plot(SF,tperiod,lag_q_BK15,'LineWidth',linew-2,'color',col_FLUX,'LineStyle',lins{iQ}); hold on;     
        end
    end
    if  (opt_admit ~= 2) % If wet or basal-flux models
        if (opt_admit == 1) 
           plot(SF,tperiod,lag_qc,'LineWidth',linew,'color',col_ECO2,'LineStyle',lins{iQ}); hold on;
        elseif (opt_admit == 3) 
           lag_qc      = data.lagBCtoSurf.qc(iQ,:);
           lag_qc_BK15 = data.lagBCtoSurf_BK15.qc(2,:);
           [~,idx_min]  = min(abs(tperiod(:)-tpmin)); 
           plot(SF,tperiod(idx_min:end),lag_qc(idx_min:end)+0.25*tperiod(idx_min:end),'LineWidth',linew,'color',col_ECO2,'LineStyle',lins{iQ}); hold on;
           plot(SF,tperiod(idx_min:end),lag_qc_BK15(idx_min:end)+0.25*tperiod(idx_min:end),'LineWidth',linew-2,'color',col_ECO2,'LineStyle',lins{iQ}); hold on;
        end
    end

    
    if (iQ==params_to_plot(1)) 
        hold('on');
        xlabel(SF,'$t_p$ [kyr]','Fontsize',fontsize,'interpreter','latex');
        set(SF,'xlim',[0 200],'Fontsize',fontsize,'Box','on');
        if (opt_admit ~= 3)
            set(SF,'ylim',[-5 20]);
            ylabel(SF,{'Lag [kyr]'},'Fontsize',fontsize,'interpreter','latex');
        else
            set(SF,'ylim',[0 100]);
            ylabel(SF,{'Bottom-to-surface lag [kyr]'},'Fontsize',fontsize,'interpreter','latex');
        end
        grid(SF,'on'); SF.XMinorGrid='on'; SF.YMinorGrid='on';
    end
     
end


%----------%%----------% Legend %----------%%----------%
switch opt_admit
    case 1 % wet
        h(1)=plot(NaN,NaN,'-','Color',col_FLUX,'linewidth',linew); hold on;
        h(2)=plot(NaN,NaN,'-','Color',col_ECO2,'linewidth',linew); hold on;
        [hh,icons,plots,txt] = legend(h,{'melt flux','carbon flux'},'Box','off','Fontsize',18,'Position',[0.28 0.82 0.2 0.1],'Units','normalized','Orientation','vertical');  
    case 2 % dry
        h(1)=plot(NaN,NaN,'-','Color',col_FLUX,'linewidth',linew); hold on;
        h(2)=plot(NaN,NaN,'-','Color',col_RMOR0,'linewidth',linew); hold on;
        h(3)=plot(NaN,NaN,'-','Color',col_RMOR1,'linewidth',linew); hold on;
        [hh,icons,plots,txt] = legend(h,{'1-d model','pseudo-2-d model ($\tau=0)$','pseudo-2-d model ($\tau>0)$'},'Box','on','EdgeColor',[1 1 1],'Fontsize',16,'Position',[0.69 0.79 0.2 0.1],'Units','normalized','Orientation','vertical');          
    case 3 % dry basal 
        h(1)=plot(NaN,NaN,'-','Color',col_ECO2,'linewidth',linew); hold on;
        h(2)=plot(NaN,NaN,'-','Color',col_ECO2,'linewidth',linew-2); hold on;
        [hh,icons,plots,txt] = legend(h,{'with reference parameters','with parameters of B\&K15'},'Box','on','EdgeColor',[1 1 1],'Fontsize',16,'Position',[0.195 0.82 0.2 0.1],'Units','normalized','Orientation','vertical');  
end
fprintf('... DONE \n\n');




