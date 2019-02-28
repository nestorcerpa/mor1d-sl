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

load('column_SLRecord.mat'); 

dt = mean(MODEL.SL.time(2:end)-MODEL.SL.time(1:end-1));
SLrecord_rate = (MODEL.SL.sl(2:end)-MODEL.SL.sl(1:end-1))/dt; % in m/kyr

%----------% Reconstructing response %----------%
gamma = MODEL.SL.dfs.gamma;
%%% phih 
MODEL.SL.dfs.gamma = gamma.*MODEL.phih_top;
COLUMN.phip_top    = idft(MODEL.SL.dfs)*MODEL.par.rhow/MODEL.par.rhom./(MODEL.par.H*1e3);
%%% clh
MODEL.SL.dfs.gamma = gamma.*MODEL.clh_top;
COLUMN.clp_top     = idft(MODEL.SL.dfs)*MODEL.par.rhow/MODEL.par.rhom./(MODEL.par.H*1e3);
%%% fh 
MODEL.SL.dfs.gamma = gamma.*MODEL.fh_top;
COLUMN.fp_top      = idft(MODEL.SL.dfs)*MODEL.par.rhow/MODEL.par.rhom./(MODEL.par.H*1e3);
%%% fch 
MODEL.SL.dfs.gamma = gamma.*MODEL.fch_top;
COLUMN.fcp_top     = idft(MODEL.SL.dfs)*MODEL.par.rhow/MODEL.par.rhom./(MODEL.par.H*1e3);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----------%%----------% FIGURES %----------%%----------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%----------%%----------% Initializing %----------%%----------%
numfig=figure();
set(numfig,'position', [0, 100, 800, 800]);

p = panel(numfig);
p.pack(6, 1);  

p.de.margin = 10;
p.margin    = [20 15 10 5];   %% margin[left bottom right top]


%----------%%----------% Options %----------%%----------%

xmin        = -800; 
xmax        = 0;
xminticks   = 17;
linew       = 3;


%----------%%----------% Plotting %----------%%----------%

%----------% S %----------%
p(1,1).select(); F11=gca;
plot(F11,MODEL.SL.time,MODEL.SL.sl,'k','linewidth',linew); hold on;
ylabel(F11,'S [m]');
set(F11,'xlim',[xmin xmax]);
set(F11,'XMinorTick','on'); F11.XAxis.MinorTickValues = linspace(xmin,xmax,xminticks);
set(F11,'XGrid','on','XMinorGrid','on','YGrid','on','YMinorGrid','off'); 

%----------% -\dot{S} %----------%
p(2,1).select(); F21=gca;
plot(F21,MODEL.SL.time(2:end),-SLrecord_rate*1e-1,'k','linewidth',linew); hold on;
ylabel(F21,'$-\dot{S}$ [cm/yr]','Interpreter','latex');
set(F21,'xlim',[xmin xmax]);
set(F21,'XMinorTick','on'); F21.XAxis.MinorTickValues = linspace(xmin,xmax,xminticks);
set(F21,'XGrid','on','XMinorGrid','on','YGrid','on','YMinorGrid','off'); 

%----------% phi %----------%
p(3,1).select(); F21=gca;
plot(F21,MODEL.SL.time,COLUMN.phip_top./MODEL.MFields.phi(end)*100.,'color',[0 1 1],'linewidth',linew); 
ylabel(F21,'$\frac{\phi^{\prime}}{\overline{\phi}}$ [$\%$]','Interpreter','latex');
set(F21,'xlim',[xmin xmax]);
set(F21,'XMinorTick','on'); F21.XAxis.MinorTickValues = linspace(xmin,xmax,xminticks);
set(F21,'XGrid','on','XMinorGrid','on','YGrid','on','YMinorGrid','off'); 

%----------% cl %----------%
p(4,1).select(); F31=gca;
plot(F31,MODEL.SL.time,COLUMN.clp_top./(MODEL.MFields.cs(end)/MODEL.par.D_vol)*100.,'color',[1 0 0],'linewidth',linew); 
ylabel(F31,'$\frac{c_l^{\prime}}{\overline{c_l}}$ [$\%$]','Interpreter','latex');
set(F31,'xlim',[xmin xmax]);
set(F31,'XMinorTick','on'); F31.XAxis.MinorTickValues = linspace(xmin,xmax,xminticks);
set(F31,'XGrid','on','XMinorGrid','on','YGrid','on','YMinorGrid','off'); 

%----------% f %----------%
p(5,1).select(); F41=gca;
plot(F41,MODEL.SL.time,COLUMN.fp_top./MODEL.MFields.f(end)*100.,'color',[0 0 1],'linewidth',linew); hold on;
ylabel(F41,'$\frac{f^{\prime}}{\overline{f}}$ [$\%$]','Interpreter','latex');
set(F41,'xlim',[xmin xmax]);
set(F41,'XMinorTick','on'); F41.XAxis.MinorTickValues = linspace(xmin,xmax,xminticks);
set(F41,'XGrid','on','XMinorGrid','on','YGrid','on','YMinorGrid','off'); 

%----------% fc %----------%
p(6,1).select(); F51=gca;
plot(F51,MODEL.SL.time,COLUMN.fcp_top./MODEL.MFields.fc(end)*100.,'color',[0.5 0 0.5],'linewidth',linew); 
xlabel(F51,'Time before present [kyr]');
ylabel(F51,'$\frac{f_c^{\prime}}{\overline{f_c}}$ [$\%$]','Interpreter','latex');
set(F51,'xlim',[xmin xmax]);
set(F51,'XMinorTick','on'); F51.XAxis.MinorTickValues = linspace(xmin,xmax,xminticks);
set(F51,'XGrid','on','XMinorGrid','on','YGrid','on','YMinorGrid','off'); 
