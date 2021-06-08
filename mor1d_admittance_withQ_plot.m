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

par = input_parameters(); % Initializing model parameters


%----------% Loading files %----------%
load('mor1d_admittance_withQ.mat');

par_array  = data.par_array;
MFieldsTop = data.MFieldsTop; 
FFieldsTop = data.FFieldsTop; 
tperiod    = data.tp_array;
nQs = size(par_array,1);

tpmin = 0.0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----------%%----------% Reading data %----------%%----------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for iQ = 1:nQs
  
    Qarray(iQ) = par_array(iQ).Q; 
    
    %----------% Extract admittance values for t_p > min kyrs %----------%
    [xmin_min,idx_min] = min(abs(tperiod(:)*1e-3-tpmin)); 
    tperiod_m     = tperiod(idx_min:end);
    RelAdm_FLUX_m = par_array(iQ,1).delta0*abs(FFieldsTop.qh(iQ,:))./MFieldsTop.q(iQ);  
    RelAdm_ECO2_m = par_array(iQ,1).delta0*abs(FFieldsTop.qch(iQ,:))./MFieldsTop.qc(iQ);
        
    %----------% Determine value and time of peak admittance %----------%
    [maxAq(iQ),idmaxAq] = max(RelAdm_FLUX_m); tmaxAq(iQ) = tperiod_m(idmaxAq);
    [maxAqc(iQ),idmaxAqc] = max(RelAdm_ECO2_m);tmaxAqc(iQ) = tperiod_m(idmaxAqc);
    
    %----------% Determine approximated value and time of peak admittance %----------%
    n        = par_array(iQ).n;
    phimax   = (par_array(iQ).Fmax/par_array(iQ).Q)^(1./n);
    w0       = par_array(iQ).W0*par_array(iQ).Fmax/phimax; % [cm/yr]
    DeltaS   = 0.1;
    Aqth(iQ)   = n*DeltaS/par_array(iQ).Hdry*par_array(iQ).rhow/par_array(iQ).rhom*w0/par_array(iQ).W0;         % Eq. (30b)
    Aqcth(iQ)  = (n-1)*DeltaS/par_array(iQ).Hdry*par_array(iQ).rhow/par_array(iQ).rhom*w0/par_array(iQ).W0;     % Eq. (30c)
    tpqth(iQ)  = n*par_array(iQ).Hdry/(w0*1e-5)*1e-3;                                                           % Eq. (23)
    tpqcth(iQ) = n*par_array(iQ).Hdry/(w0*1e-5)*1e-3;                                                           % Eq. (23)
        
end

%----------% Extracting prefactor for fitting data at \mathcal{Q} = 10^5. %----------%
[~,idx_Q1e5] = min(abs(Qarray-1e5));
Bq  = (2*maxAq(idx_Q1e5))/Aqth(idx_Q1e5);
Bqc = (2*maxAqc(idx_Q1e5))/Aqcth(idx_Q1e5);
Cq  = (tmaxAq(idx_Q1e5))/tpqth(idx_Q1e5);
Cqc = (tmaxAqc(idx_Q1e5))/tpqcth(idx_Q1e5);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----------%%----------% FIGURE %----------%%----------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nfig = 0; 

%----------%%----------% Options %----------%%----------%

linew = 3;      % Line width for plots
legendsize= 16; % Legend fontsize
fontsize = 16;  % Labels fontsize
marksize = 12;

col_FLUX = [0 0 1];
col_ECO2 = [0.5 0 0.5];

%----------%%----------% Initializing %----------%%----------% 
nfig = nfig + 1; figure(nfig); set(nfig, 'Position', [500, 100, 650, 350]); fig1=nfig;

p = panel(nfig);

% Create panel 2x3
p.pack(2,1);

p.margin=[22 20 5 8]; %% margin[left bottom right top]
p.de.margin = 10;

p = panel(nfig);

% Create panel
p.pack(1, {99/100 1/100});
p(1,1).pack(1,2)


p.margin=[16 20 2 5]; %% margin[left bottom right top]
p.de.margin = 0;
p(1,1,1,1).marginright = 24;


%----------%%----------% Plotting max(A) %----------%%----------% 
Amaxlims = [0 45];

p(1,1,1,1).select(); SF=gca; 

%----------% Plot calculated values %----------%
semilogx(SF,Qarray,2.*maxAq*100.,'o','markersize',marksize,'color',col_FLUX,'MarkerFaceColor',col_FLUX); hold on; 
semilogx(SF,Qarray,2.*maxAqc*100.,'^','markersize',marksize,'color',col_ECO2,'MarkerFaceColor',col_ECO2); hold on; 
%----------% Plot theoretical values %----------%
semilogx(SF,Qarray,Bq*Aqth*100,'-','linewidth',linew,'color',col_FLUX); hold on;   leg1{1} = sprintf('$B_{\\mathcal{Q}}$=%4.2f',Bq);
semilogx(SF,Qarray,Bqc*Aqcth*100,'-','linewidth',linew,'color',col_ECO2); hold on; leg1{2} = sprintf('$B_{\\mathcal{Q}_c}$=%4.2f',Bqc);

set(SF,'Fontsize',fontsize,'Box','on');
xlabel(SF,'$\mathcal{Q}$')
ylabel(SF,'max($A$) [$\%$]');
set(SF,'xlim',[1e4 1e6]);
set(SF,'YTick',[Amaxlims(1):10:Amaxlims(2)]); SF.YAxis.MinorTick='on'; SF.YAxis.MinorTickValues=Amaxlims(1):5:Amaxlims(2);
grid(SF,'on'); SF.XMinorGrid='on'; SF.YMinorGrid='on';
text(1e5,28,leg1{1},'fontsize',fontsize+1);
text(2.5e5,12,leg1{2},'fontsize',fontsize+1);


%----------%%----------% Plotting time(max(A)) %----------%%----------% 
tmaxlims = [0 160];

p(1,1,1,2).select(); SF=gca; 

%----------% Plot calculated values %----------%
semilogx(SF,Qarray,tmaxAq,'o','markersize',marksize,'color',col_FLUX,'MarkerFaceColor',col_FLUX); hold on; 
semilogx(SF,Qarray,tmaxAqc,'^','markersize',marksize,'color',col_ECO2,'MarkerFaceColor',col_ECO2); hold on; 
%----------% Plot theoretical values %----------%
semilogx(SF,Qarray,Cq*tpqth,'-','linewidth',linew,'color',col_FLUX); hold on;      leg1{1} = sprintf('$C_{\\mathcal{Q}}$=%4.2f',Cq);
semilogx(SF,Qarray,Cqc*tpqcth,'-','linewidth',linew,'color',col_ECO2); hold on;    leg1{2} = sprintf('$C_{\\mathcal{Q}_c}$=%4.2f',Cqc);

set(SF,'Fontsize',fontsize,'Box','on');
xlabel(SF,'$\mathcal{Q}$')
ylabel(SF,'period$({\textrm{max}(A)})$ [kyr]');
set(SF,'xlim',[1e4 1e6]);
set(SF,'YTick',[tmaxlims(1):20:tmaxlims(2)]); SF.YAxis.MinorTick='on'; SF.YAxis.MinorTickValues=tmaxlims(1):10:tmaxlims(2);
grid(SF,'on'); SF.XMinorGrid='on'; SF.YMinorGrid='on';
text(1.02e4,50,leg1{1},'fontsize',fontsize+1);
text(2.2e4,85,leg1{2},'fontsize',fontsize+1);


%----------%%----------% Legend %----------%%----------%
p(1,2).select(); F01=gca; F01.Visible = 'off'; hp=zeros(3,1);
h(1)=plot(NaN,NaN,'-','Color',col_FLUX,'linewidth',linew); hold on;
h(2)=plot(NaN,NaN,'-','Color',col_ECO2,'linewidth',linew); hold on;
[hh,icons,plots,txt] = legend(h,{'melt flux','carbon flux'},'Box','off','Fontsize',18,'Position',[0.08 0.82 0.2 0.1],'Units','normalized','Orientation','vertical');  

