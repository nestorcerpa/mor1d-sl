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

clear all; close all;

%----------% Add paths to functions %----------%
restoredefaultpath;
addpath(genpath([pwd,'/src/']),genpath([pwd,'/external-functions/'])); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----------%   Loading data files     %----------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('mor1d_SLRecord.mat'); 

dt = mean(DATA.SL.time(2:end)-DATA.SL.time(1:end-1));
SLrecord_rate = (DATA.SL.sl(2:end)-DATA.SL.sl(1:end-1))/dt; % in m/kyr

%%% Reorder models in increasing order of w_0
iw_w0ref = 1; % Results with reference w0 is computed first in mor1d_slrecord_compute.m; 
for iw=1:length(DATA.model) 
    Q_values(iw) = DATA.model{iw}.par.Q;
    TMP{iw} = DATA.model{iw};
end
[~,idx_sort] = sort(Q_values); iw_w0ref = idx_sort(iw_w0ref);
for iw=1:length(DATA.model)
   DATA.model{iw} = TMP{idx_sort(iw)};
end;


%----------% Reconstructing response %----------%
fprintf('\nReconstructing SL-response...');
gamma = DATA.SL.dfs.gamma;
for iw=1:length(DATA.model)
    fprintf('\n---> model : %2d',iw);
    %%% phih 
    DATA.SL.dfs.gamma  = gamma.*DATA.model{iw}.phih_top;
    COLUMN{iw}.phip_top= idft(DATA.SL.dfs)*DATA.model{iw}.par.rhow/DATA.model{iw}.par.rhom./(DATA.model{iw}.par.H*1e3);
    %%% clh
    DATA.SL.dfs.gamma  = gamma.*DATA.model{iw}.clh_top;
    COLUMN{iw}.clp_top = idft(DATA.SL.dfs)*DATA.model{iw}.par.rhow/DATA.model{iw}.par.rhom./(DATA.model{iw}.par.H*1e3);
    %%% fh 
    DATA.SL.dfs.gamma  = gamma.*DATA.model{iw}.qh_top;
    COLUMN{iw}.qp_top  = idft(DATA.SL.dfs)*DATA.model{iw}.par.rhow/DATA.model{iw}.par.rhom./(DATA.model{iw}.par.H*1e3);
    %%% fch 
    DATA.SL.dfs.gamma  = gamma.*DATA.model{iw}.qch_top;
    COLUMN{iw}.qcp_top = idft(DATA.SL.dfs)*DATA.model{iw}.par.rhow/DATA.model{iw}.par.rhom./(DATA.model{iw}.par.H*1e3);
    
    %%% Rmorh0 
    DATA.SL.dfs.gamma  = gamma.*DATA.model{iw}.Rmorh0;
    COLUMN{iw}.Rmor0   = idft(DATA.SL.dfs)*DATA.model{iw}.par.rhow/DATA.model{iw}.par.rhom./(DATA.model{iw}.par.H*1e3);
    %%% Rmorh1 
    DATA.SL.dfs.gamma  = gamma.*DATA.model{iw}.Rmorh1;
    COLUMN{iw}.Rmor1   = idft(DATA.SL.dfs)*DATA.model{iw}.par.rhow/DATA.model{iw}.par.rhom./(DATA.model{iw}.par.H*1e3);
    
    %%% Rcmorh0 
    DATA.SL.dfs.gamma  = gamma.*DATA.model{iw}.Rcmorh0;
    COLUMN{iw}.Rcmor0  = idft(DATA.SL.dfs)*DATA.model{iw}.par.rhow/DATA.model{iw}.par.rhom./(DATA.model{iw}.par.H*1e3);
    %%% Rcmorh1 
    DATA.SL.dfs.gamma  = gamma.*DATA.model{iw}.Rcmorh1;
    COLUMN{iw}.Rcmor1  = idft(DATA.SL.dfs)*DATA.model{iw}.par.rhow/DATA.model{iw}.par.rhom./(DATA.model{iw}.par.H*1e3);    
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----------%%----------% FIGURE Time series %----------%%----------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

col_RMOR0  = [0.4 0.6 1.0];
col_RMOR1  = [0.0 0.2 0.4];
col_RCMOR0 = [1.0 0.6 1.0];
col_RCMOR1 = [0.8 0.0 0.4];

iw_to_plot = iw_w0ref; par = DATA.model{iw_to_plot}.par; w0_ref = par.Fmax^((par.n-1)/par.n) * par.Q^(1./par.n)*par.W0; 
fprintf('\n -----> Plotting results for case max(w) = %4.1f m/yr \n',w0_ref*0.01)


%----------%%----------% Initializing %----------%%----------%
numfig=figure();
set(numfig,'position', [0, 100, 800, 800]);

p = panel(numfig);
p.pack(5, 1);  

p.fontsize  = 16;

p.de.margin = 5;
p.margin    = [22 16 10 5];   %% margin[left bottom right top]


%----------%%----------% Options %----------%%----------%

xmin        = -800;  [~,idx_min] = min(abs(DATA.SL.time-xmin));
xmax        = 0;     [~,idx_max] = min(abs(DATA.SL.time-xmax));
xminticks   = 17;  
linew       = 3;


%----------%%----------% Plotting %----------%%----------%

%----------% S %----------%
p(1,1).select(); F11=gca;
plot(F11,DATA.SL.time,DATA.SL.sl,'k','linewidth',linew); hold on;
ylabel(F11,'$S$ [m]','Interpreter','latex');
set(F11,'xlim',[xmin xmax],'XtickLabel',[]);
set(F11,'XMinorTick','on'); F11.XAxis.MinorTickValues = linspace(xmin,xmax,xminticks);
set(F11,'ylim',[-75,75]);
set(F11,'YMinorTick','on'); F11.YAxis.MinorTickValues = linspace(-50,50,5);
set(F11,'XGrid','on','XMinorGrid','on','YGrid','on','YMinorGrid','on'); set(F11,'Box','on');


%----------% -\dot{S} %----------%
p(2,1).select(); F21=gca;
plot(F21,DATA.SL.time(2:end),-SLrecord_rate*1e-1,'k','linewidth',linew); hold on;
ylabel(F21,'$-\dot{S}$ [cm/yr]','Interpreter','latex');
set(F21,'xlim',[xmin xmax],'XtickLabel',[]);
set(F21,'XMinorTick','on'); F21.XAxis.MinorTickValues = linspace(xmin,xmax,xminticks);
set(F21,'ylim',[-1.5,1.5]);
set(F21,'YMinorTick','on'); F21.YAxis.MinorTickValues = linspace(-1.5,1.5,7);
set(F21,'XGrid','on','XMinorGrid','on','YGrid','on','YMinorGrid','on'); set(F21,'Box','on');


%----------% f %----------%
p(3,1).select(); F51=gca;
plot(F51,DATA.SL.time,COLUMN{iw_to_plot}.Rmor0./DATA.model{iw_to_plot}.MFields.Rmor*100.,'-','color',col_RMOR0,'linewidth',linew); hold on;
plot(F51,DATA.SL.time,COLUMN{iw_to_plot}.qp_top./DATA.model{iw_to_plot}.MFields.q(end)*100.,'color',[0 0 1],'linewidth',linew); hold on;
plot(F51,DATA.SL.time,COLUMN{iw_to_plot}.Rmor1./DATA.model{iw_to_plot}.MFields.Rmor*100.,'-','color',col_RMOR1,'linewidth',linew); hold on;
ylabel(F51,'$\frac{|Q^{\prime}|}{\overline{Q}}$ [$\%$]','Interpreter','latex');
set(F51,'xlim',[xmin xmax],'XtickLabel',[]);
set(F51,'XMinorTick','on'); F51.XAxis.MinorTickValues = linspace(xmin,xmax,xminticks);
set(F51,'ylim',[-10,10]);
set(F51,'YMinorTick','on'); F51.YAxis.MinorTickValues = linspace(-10,10,5);
set(F51,'XGrid','on','XMinorGrid','on','YGrid','on','YMinorGrid','on'); set(F51,'Box','on');


%----------% fc %----------%
p(4,1).select(); F61=gca;
plot(F61,DATA.SL.time,COLUMN{iw_to_plot}.Rcmor0./DATA.model{iw_to_plot}.MFields.Rcmor*100.,'-','color',col_RCMOR0,'linewidth',linew); hold on;
plot(F61,DATA.SL.time,COLUMN{iw_to_plot}.qcp_top./DATA.model{iw_to_plot}.MFields.qc(end)*100.,'color',[0.5 0 0.5],'linewidth',linew); hold on;
plot(F61,DATA.SL.time,COLUMN{iw_to_plot}.Rcmor1./DATA.model{iw_to_plot}.MFields.Rcmor*100.,'-','color',col_RCMOR1,'linewidth',linew); 
xlabel(F61,'Time before present [kyr]');
ylabel(F61,'$\frac{|Q_c^{\prime}|}{\overline{Q}_c}$ [$\%$]','Interpreter','latex');
set(F61,'xlim',[xmin xmax]);
set(F61,'XMinorTick','on'); F61.XAxis.MinorTickValues = linspace(xmin,xmax,xminticks);
set(F61,'ylim',[-7,7]);
set(F61,'YMinorTick','on'); F61.YAxis.MinorTickValues = linspace(-5,5,5);
set(F61,'XGrid','on','XMinorGrid','on','YGrid','on','YMinorGrid','on'); set(F61,'Box','on');




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----------%%----------% FIGURE Cross-Correlation %----------%%----------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%----------% Cross-correlation q %----------%
[xc_q,lags_q] = xcorr(COLUMN{iw_to_plot}.qp_top(2:end),-SLrecord_rate);
[~,I] = max(abs(xc_q)); fprintf('\n -> Peak lag q at %4.2f kyr \n',lags_q(I)*dt);
 
%----------% Cross-correlation qc %----------%
[xc_qc,lags_qc] = xcorr(COLUMN{iw_to_plot}.qcp_top(2:end),-SLrecord_rate);
[~,I] = max(abs(xc_qc)); fprintf(' -> Peak lag qc at %4.2f kyr \n',lags_qc(I)*dt);

%----------% Cross-correlation Rmor0 %----------%
[xc_Rmor0,lags_Rmor0] = xcorr(COLUMN{iw_to_plot}.Rmor0(2:end),-SLrecord_rate);
[~,I] = max(abs(xc_Rmor0)); fprintf(' -> Peak lag Rmor0 at %4.2f kyr \n',lags_Rmor0(I)*dt);

%----------% Cross-correlation Rmor1 %----------%
[xc_Rmor1,lags_Rmor1] = xcorr(COLUMN{iw_to_plot}.Rmor1(2:end),-SLrecord_rate);
[~,I] = max(abs(xc_Rmor1)); fprintf(' -> Peak lag Rmor1 at %4.2f kyr \n',lags_Rmor1(I)*dt);

%----------% Cross-correlation Rcmor0 %----------%
[xc_Rcmor0,lags_Rcmor0] = xcorr(COLUMN{iw_to_plot}.Rcmor0(2:end),-SLrecord_rate);
[~,I] = max(abs(xc_Rcmor0)); fprintf(' -> Peak lag Rcmor0 at %4.2f kyr \n',lags_Rcmor0(I)*dt);

%----------% Cross-correlation Rcmor1 %----------%
[xc_Rcmor1,lags_Rcmor1] = xcorr(COLUMN{iw_to_plot}.Rcmor1(2:end),-SLrecord_rate);
[~,I] = max(abs(xc_Rcmor1)); fprintf(' -> Peak lag Rcmor1 at %4.2f kyr \n\n',lags_Rcmor1(I)*dt);


%----------%%----------% Options %----------%%----------%
xmin_large  = 0; 
xmax_large  = 140;

p(5,1).pack(1, 2);  
p(5).margintop = 18;
p(5,1,1,1).marginright = 10;
p(5,1,1,2).marginleft = 10;

%----------%%----------% Plot Cross-correlation for Q %----------%%----------%
p(5,1,1,1).select(); F00=gca;
hf{1}=plot(F00,lags_Rmor0*dt,xc_Rmor0/max(abs(xc_Rmor0)),'color',col_RMOR0,'linewidth',linew); hold on;
hf{2}=plot(F00,lags_q*dt,xc_q/max(abs(xc_q)),'color',[0 0 1],'linewidth',linew); hold on;
hf{3}=plot(F00,lags_Rmor1*dt,xc_Rmor1/max(abs(xc_Rmor1)),'color',col_RMOR1,'linewidth',linew); hold on;
xlabel(F00,'Lags [kyr]');
ylabel(F00,{'Cross-correlation'},'Interpreter','latex');
set(F00,'xlim',[xmin_large xmax_large],'Xtick',[0:20:140]);
set(F00,'XMinorTick','on'); F00.XAxis.MinorTickValues = linspace(xmin_large,xmax_large,29);
set(F00,'ylim',[-1.0,1.0]);
set(F00,'YMinorTick','on'); F00.YAxis.MinorTickValues = linspace(-1.0,1.0,21);
set(F00,'XGrid','on','XMinorGrid','on','YGrid','on','YMinorGrid','on'); set(F00,'Box','on');
h = legend([hf{2},hf{1},hf{3}],{'1-d model','pseudo-2-d model ($\tau=0$)','pseudo-2-d model ($\tau>0$)',},'location','southeast','fontsize',12);
text(F00,50,0.85,'Melt flux','fontsize',18);

%----------%%----------% Plot Cross-correlation for Q %----------%%----------%
p(5,1,1,2).select(); F00=gca;
hf{1}=plot(F00,lags_Rcmor0*dt,xc_Rcmor0/max(abs(xc_Rcmor0)),'color',col_RCMOR0,'linewidth',linew); hold on;
hf{2}=plot(F00,lags_qc*dt,xc_qc/max(abs(xc_qc)),'color',[.5 0 .5],'linewidth',linew); hold on;
hf{3}=plot(F00,lags_Rcmor1*dt,xc_Rcmor1/max(abs(xc_Rcmor1)),'color',col_RCMOR1,'linewidth',linew); hold on;
xlabel(F00,'Lags [kyr]');
%ylabel(F00,{'Cross-correlation'},'Interpreter','latex');
set(F00,'xlim',[xmin_large xmax_large],'Xtick',[0:20:140]);
set(F00,'XMinorTick','on'); F00.XAxis.MinorTickValues = linspace(xmin_large,xmax_large,29);
set(F00,'ylim',[-1.0,1.0]);
set(F00,'YMinorTick','on'); F00.YAxis.MinorTickValues = linspace(-1.0,1.0,21);
set(F00,'XGrid','on','XMinorGrid','on','YGrid','on','YMinorGrid','on'); set(F00,'Box','on');
set(F00,'YAxisLocation','right');
h = legend([hf{2},hf{1},hf{3}],{'1-d model','pseudo-2-d model ($\tau=0$)','pseudo-2-d model ($\tau>0$)',},'location','southeast','fontsize',12);
text(F00,50,0.85,'Carbon flux','fontsize',18);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----------%%----------% FIGURE Sensitivity to w_0 %----------%%----------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%----------%%----------% Initializing %----------%%----------%
numfig=figure(); 
set(numfig,'position', [0, 100, 700, 500]);

p = panel(numfig);
p.pack(2, 2);  

p.fontsize  = 16;

p.de.margin = 8;
p.margin    = [16 16 5 5];   %% margin[left bottom right top]
p(1,1).marginright = 25;
p(2,1).marginright = 25;

%----------%%----------% Melt Flux %----------%%----------%

p(1,1).select(); fig=gca;
hold(fig,'on');
for iw = 1:length(DATA.model)
    par = DATA.model{iw}.par; 
    w0(iw) = par.Fmax^((par.n-1)/par.n) * par.Q^(1./par.n)*par.W0;
    Q   = COLUMN{iw}.qp_top./DATA.model{iw}.MFields.q(end)*100.;
    R0  = COLUMN{iw}.Rmor0./DATA.model{iw}.MFields.Rmor(end)*100.;
    R1  = COLUMN{iw}.Rmor1./DATA.model{iw}.MFields.Rmor(end)*100.;
    
    varQ(iw)  = max(Q(idx_min:idx_max))-min(Q(idx_min:idx_max));
    varR0(iw) = max(R0(idx_min:idx_max))-min(R0(idx_min:idx_max));
    varR1(iw) = max(R1(idx_min:idx_max))-min(R1(idx_min:idx_max));
    
    if iw==iw_w0ref
       w0ref=w0(iw); varQref = varQ(iw); varR0ref = varR0(iw); varR1ref = varR1(iw);
    end
end
plot(w0*1e-2,varR0,'-','color',col_RMOR0,'linewidth',linew);
plot(w0*1e-2,varQ,'-','color','b','linewidth',linew); 
plot(w0*1e-2,varR1,'-','color',col_RMOR1,'linewidth',linew);
plot(w0ref*1e-2,varR0ref,'o','markersize',10,'MarkerEdgeColor',col_RMOR0,'MarkerFaceColor',col_RMOR0); 
plot(w0ref*1e-2,varQref,'o','markersize',10,'MarkerEdgeColor','b','MarkerFaceColor','b');
plot(w0ref*1e-2,varR1ref,'o','markersize',10,'MarkerEdgeColor',col_RMOR1,'MarkerFaceColor',col_RMOR1); 
xlabel(fig,'Maximum melt velocity [m/yr]');
ylabel(fig,'$A_Q$ [$\%$]');
xlim(fig,[0 10]);
ylim(fig,[0 40]);
set(fig,'Box','on','XtickLabel',[]);
grid(fig,'on');


p(1,2).select(); fig=gca;
hold(fig,'on');
for iw = 1:length(DATA.model)
    
    par = DATA.model{iw}.par; 
    w0(iw) = par.Fmax^((par.n-1)/par.n) * par.Q^(1./par.n)*par.W0;
    %fprintf('w0=%5.1f [m/y]',w0*1e-2)
    %----------% Cross-correlation f %----------%
    [xc_q,lags] = xcorr(COLUMN{iw}.qp_top(2:end),-SLrecord_rate);
    [~,I] = max(abs(xc_q)); 
    lag_q(iw) = lags(I);
    %----------% Cross-correlation Rmor0 %----------%
    [xc_Rmor0,lags] = xcorr(COLUMN{iw}.Rmor0(2:end),-SLrecord_rate);
    [~,IR0] = max(abs(xc_Rmor0));
    lag_Rmor0(iw) = lags(IR0);
    %----------% Cross-correlation Rmor0 %----------%
    [xc_Rmor1,lags] = xcorr(COLUMN{iw}.Rmor1(2:end),-SLrecord_rate);
    [~,IR1] = max(abs(xc_Rmor1));
    lag_Rmor1(iw) = lags(IR1);
    
    if iw==iw_w0ref
        w0ref=w0(iw); lag_qref = lag_q(iw); lag_Rmor0ref=lag_Rmor0(iw); lag_Rmor1ref=lag_Rmor1(iw); 
    end
end
plot(w0*1e-2,lag_Rmor0*dt,'-','Color',col_RMOR0,'linewidth',linew);
plot(w0*1e-2,lag_q*dt,'-','Color','b','linewidth',linew); 
plot(w0*1e-2,lag_Rmor1*dt,'-','Color',col_RMOR1,'linewidth',linew);
plot(w0ref*1e-2,lag_Rmor0ref*dt,'o','markersize',10,'MarkerEdgeColor',col_RMOR0,'MarkerFaceColor',col_RMOR0); 
plot(w0ref*1e-2,lag_qref*dt,'o','markersize',10,'MarkerEdgeColor','b','MarkerFaceColor','b'); 
plot(w0ref*1e-2,lag_Rmor1ref*dt,'o','markersize',10,'MarkerEdgeColor',col_RMOR1,'MarkerFaceColor',col_RMOR1); 
xlabel(fig,'Maximum melt velocity [m/yr]');
ylabel(fig,'Lag [kyr]');
xlim(fig,[0 10]);
ylim(fig,[0 25]);
set(fig,'Box','on','XtickLabel',[]);
grid(fig,'on');



%----------%%----------% Carbon Flux %----------%%----------%

p(2,1).select(); fig=gca;
hold(fig,'on');
for iw = 1:length(DATA.model)
    
    par = DATA.model{iw}.par; 
    w0(iw) = par.Fmax^((par.n-1)/par.n) * par.Q^(1./par.n)*par.W0;
    
    Qc = COLUMN{iw}.qcp_top./DATA.model{iw}.MFields.qc(end)*100.;
    Rc0 = COLUMN{iw}.Rcmor0./DATA.model{iw}.MFields.Rcmor(end)*100.;
    Rc1 = COLUMN{iw}.Rcmor1./DATA.model{iw}.MFields.Rcmor(end)*100.;
    
    varQc(iw)   = max(Qc(idx_min:idx_max))-min(Qc(idx_min:idx_max));
    varRc0(iw) = max(Rc0(idx_min:idx_max))-min(Rc0(idx_min:idx_max));
    varRc1(iw) = max(Rc1(idx_min:idx_max))-min(Rc1(idx_min:idx_max));
    
    if iw==iw_w0ref
       w0ref=w0(iw); varQcref = varQc(iw); varRc0ref = varRc0(iw); varRc1ref = varRc1(iw);
    end
    
end
plot(w0*1e-2,varRc0,'-','Color',col_RCMOR0,'linewidth',linew);
plot(w0*1e-2,varQc,'-','Color','m','linewidth',linew);
plot(w0*1e-2,varRc1,'-','Color',col_RCMOR1,'linewidth',linew);       
plot(w0ref*1e-2,varRc0ref,'o','markersize',10,'MarkerEdgeColor',col_RCMOR0,'MarkerFaceColor',col_RCMOR0);
plot(w0ref*1e-2,varQcref,'o','markersize',10,'MarkerEdgeColor','m','MarkerFaceColor','m');
plot(w0ref*1e-2,varRc1ref,'o','markersize',10,'MarkerEdgeColor',col_RCMOR1,'MarkerFaceColor',col_RCMOR1);
xlabel(fig,'Maximum melt velocity [m/yr]');
ylabel(fig,'$A_{Q_c}$ [$\%$]');
xlim(fig,[0 10]);
ylim(fig,[0 40]);
set(fig,'Box','on');
grid(fig,'on');



p(2,2).select(); fig=gca;
hold(fig,'on');
for iw = 1:length(DATA.model)
    
    par = DATA.model{iw}.par; 
    w0(iw) = par.Fmax^((par.n-1)/par.n) * par.Q^(1./par.n)*par.W0;
    %fprintf('w0=%5.1f [m/y]',w0*1e-2)

    %----------% Cross-correlation qc %----------%
    [xc_q,lags] = xcorr(COLUMN{iw}.qcp_top(2:end),-SLrecord_rate);
    [~,I] = max(abs(xc_q)); 
    lag_qc(iw) = lags(I);
    %----------% Cross-correlation Rcmor0 %----------%
    [xc_Rmor0,lags] = xcorr(COLUMN{iw}.Rcmor0(2:end),-SLrecord_rate);
    [~,IR0] = max(abs(xc_Rmor0)); 
    lag_Rcmor0(iw) = lags(IR0);
    %----------% Cross-correlation Rcmor0 %----------%
    [xc_Rmor1,lags] = xcorr(COLUMN{iw}.Rcmor1(2:end),-SLrecord_rate);
    [~,IR1] = max(abs(xc_Rmor1));
    lag_Rcmor1(iw) = lags(IR1);

    if iw==iw_w0ref
        w0ref=w0(iw); lag_qcref = lag_qc(iw); lag_Rcmor0ref=lag_Rcmor0(iw); lag_Rcmor1ref=lag_Rcmor1(iw); 
    end
    
end
plot(w0*1e-2,lag_Rcmor0*dt,'-','Color',col_RCMOR0,'linewidth',linew);
plot(w0*1e-2,lag_qc*dt,'-','Color','m','linewidth',linew); 
plot(w0*1e-2,lag_Rcmor1*dt,'-','Color',col_RCMOR1,'linewidth',linew);
plot(w0ref*1e-2,lag_Rcmor0ref*dt,'o','markersize',10,'MarkerEdgeColor',col_RCMOR0,'MarkerFaceColor',col_RCMOR0); 
plot(w0ref*1e-2,lag_qcref*dt,'o','markersize',10,'MarkerEdgeColor','m','MarkerFaceColor','m'); 
plot(w0ref*1e-2,lag_Rcmor1ref*dt,'o','markersize',10,'MarkerEdgeColor',col_RCMOR1,'MarkerFaceColor',col_RCMOR1); 
xlabel(fig,'Maximum melt velocity [m/yr]');
ylabel(fig,'Lag [kyr]');
xlim(fig,[0 10]);
ylim(fig,[0 25]);
set(fig,'Box','on');
grid(fig,'on');

