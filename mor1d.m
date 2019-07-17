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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----------% Initialize parameters %----------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

par = input_parameters(); % Initializing model parameters

%par.H = par.Hdry; 
%par   = get_dimensionless_parameters(par);  % Updating dimensionless parameters 

%----------% Spatial and time arrays %----------%
zarray = linspace(0,1,par.nz);

%----------% Plotting parameters %----------%
linew = 3;  % linewidth
fonts = 18; % fontsize

nperiods = 2; % number of fluctuating-periods to plot in z-t space

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----------%    Solving problem     %----------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%----------% Calculate steady state %----------%
[MFields.cs,MFields.phi,~] = mean_analytical(zarray,par);
%----------% Get other steady-state variables %----------%
[MFields.W,MFields.q,MFields.qc] = get_other_mfields(MFields.phi,MFields.cs,par);
    
%----------% Calculate fluctuations %----------%
[FFields.csh,FFields.phih,~] = fluctuations(zarray,par);
%----------% Get other fluctuating variables %----------%
[FFields.Wh,FFields.qh,FFields.qch] = get_other_ffields(MFields.phi,MFields.cs,FFields.phih,FFields.csh,par); 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----------%  2-d Approximation for MOR focusing  %----------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

drho = 500;           % solid-liquid density contrast 
rhol = par.rhom-drho; % liquid density
rhoc = 2900.;         % crust density

par.alpha = 30.;      % angle of decompaction channel
par.xf    = 33.5;     % focusing distance xf
hf    = par.xf*tan(par.alpha*pi/180.); % depth decompating channel at x=xf

nz = 30;
zf_array  = linspace(par.H-hf,par.H,nz)/par.H;
tan_alpha = tan(par.alpha*pi/180.);

%----------% Calculate steady state %----------%
MFields.Rmor = 0.0; MFields.Rcmor = 0.0; MFields.Hcru = 0.0; tau_array(1) = nan;
for i=2:length(zf_array)
    dz  = zf_array(i)-zf_array(i-1);
    [cs_i,phi_i,~] = mean_analytical(zf_array(i),par);  % Calculate base state analytically 
    [W_i,q_i,qc_i] = get_other_mfields(phi_i,cs_i,par); % Get other variables    
    MFields.Rmor  = MFields.Rmor  + 2*q_i*dz/tan_alpha; % Factor 2 is for taking into account both sides of ridge's axis
    MFields.Rcmor = MFields.Rcmor + 2*qc_i*dz/tan_alpha;
    MFields.Hcru  = MFields.Hcru  + 2*q_i/(pi/2)*rhol/rhoc*dz/tan_alpha;
    tau_array(i)  = 0.0; % tranport-time for melt originated at z = zf_array(i)
    for j = i:length(zf_array)
        dzj  = zf_array(j)-zf_array(j-1);
        [cs_j,phi_j,~] = mean_analytical(zf_array(j),par);  % Calculate base state analytically 
        [W_j,f_j,fc_j] = get_other_mfields(phi_j,cs_j,par); % Get other variables   
        w_j = par.Q*phi_j^(par.n-1)*(1-phi_j)^2;            % steady-state melt velocity
        tau_array(i) = tau_array(i) + 1./w_j*dzj;
    end
end
fprintf("\n ---> Mean crustal thickness = %5.2f km \n",MFields.Hcru*par.H);


%----------% Calculate fluctuations %----------%
par.verb = 'off';
fprintf('\n Calculating approx 2-d model (complex) ...\n');
FFields.Rmorh0 = 0.0; FFields.Rcmorh0 = 0.0;
FFields.Rmorh1 = 0.0; FFields.Rcmorh1 = 0.0;
for i=2:length(zf_array)
    dz  = zf_array(i)-zf_array(i-1);
    [cs_i,phi_i,~] = mean_analytical(zf_array(i),par); % Calculate base state analytically 
    [ch_i,phih_i] = fluctuations(zf_array(i),par);
    [Wh_i,qh_i,qch_i] = get_other_ffields(phi_i,cs_i,phih_i,ch_i,par); 
    % Instantaneous focusing     
    FFields.Rmorh0  = FFields.Rmorh0  + 2.0*qh_i*dz/tan_alpha;
    FFields.Rcmorh0 = FFields.Rcmorh0 + 2.0*qch_i*dz/tan_alpha;
    % Focusing time equal to steady-state melt transport time
    %w_i = par.Q*phi_i^(par.n-1)*(1-phi_i)^2;           % steady-state melt velocity
    tau =  tau_array(i); %(1-zf_array(i))/tan_alpha/w_i;
    FFields.Rmorh1   = FFields.Rmorh1  + 2.0*exp(-1i*par.omega*tau)*qh_i*dz/tan_alpha;
    FFields.Rcmorh1  = FFields.Rcmorh1 + 2.0*exp(-1i*par.omega*tau)*qch_i*dz/tan_alpha;
end
fprintf('\n... Done \n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----------%   Admittance calculation   %----------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

AQ  = 2.*par.delta0*abs(FFields.qh(end))/MFields.q(end); 
AR0 = 2.*par.delta0*abs(FFields.Rmorh0)/MFields.Rmor; 
AR1 = 2.*par.delta0*abs(FFields.Rmorh1)/MFields.Rmor;

fprintf('\n \t ..Admittance flux per 100-m of SL change in 1-d model           : %5.2f %%',AQ*100.)
fprintf('\n \t ..Admittance flux per 100-m of SL change in 2-d approx. (tau=0) : %5.2f %%',AR0*100.)
fprintf('\n \t ..Admittance flux per 100-m of SL change in 2-d approx. (tau>0) : %5.2f %% \n',AR1*100.)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----------%    Plotting solution   %----------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nfig  = 0; 

%----------% Plot mean state %----------%
nfig=nfig+1; figure(nfig); 
plot_meanfields(nfig,MFields,zarray,par,linew,fonts); 

%----------% Plot fluctuations %----------%
nfig=nfig+1; figure(nfig); 
plot_fluctuationsh(nfig,FFields,zarray,par,linew,fonts);

nfig=nfig+1; figure(nfig);
plot_fluctuations_ztspace(nfig,MFields,FFields,par,linew,fonts,nperiods); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%         LOCAL FUNCTIONS         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_meanfields(nfig,MFields,zarray,par,linew,fontsize)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT_MEANFIELDS Plots mean variables \bar{c},\bar{phi},\bar{f} and \bar{f}_c
%
%   n.b.: this function uses external-package panel-2.12
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %----------%    Figure  Properties    %----------%    
    set(nfig,'position',[0 50 600 300]); 
    p=panel(nfig); p.pack(1,4); 
    p.fontsize=fontsize;
    p.title('Mean state')
       
    p.margin    = [20 20 15 10]; % Margins outside [left bottom right top]
    p.de.margin = 11;            % Margins between sub-panels
    
    
    %----------% Porosity %----------%
    p(1,1).select(); F01=gca; 
    semilogx(F01,MFields.phi,zarray,'k','linewidth',linew); hold(F01,'on');
    xlabel(F01,'$\overline{\phi}$','Fontsize',fontsize);
    ylabel(F01,'$z$'); 
    set(F01,'xlim',[1e-8 1e-2],'Xtick',[1e-7 1e-5 1e-3])
    set(F01,'ylim',[0 1],'Ytick',0:0.25:1)
    set(F01,'yaxislocation','left');       
    set(F01,'Box','on');
    grid(F01,'on'); set(F01,'XMinorGrid','on','YMinorGrid','on'); 
      
    %----------% Concentration in solid %----------%
    p(1,2).select(); F01=gca; 
    semilogx(F01,MFields.cs,zarray,'k','linewidth',linew); hold(F01,'on');
    xlabel(F01,'$\overline{c}$','Fontsize',fontsize);
    set(F01,'ylim',[0 1],'Ytick',0:0.25:1)
    set(F01,'YtickLabel',[]);       
    set(F01,'Box','on');
    grid(F01,'on'); set(F01,'XMinorGrid','on','YMinorGrid','on'); 
    
    %----------% Melt flux %----------%
    p(1,3).select(); F01=gca; 
    semilogx(F01,MFields.q,zarray,'k','linewidth',linew); hold(F01,'on');
    xlabel(F01,'$\overline{Q}$','Fontsize',fontsize);
    set(F01,'ylim',[0 1],'Ytick',0:0.25:1)
    set(F01,'YtickLabel',[]);       
    set(F01,'Box','on');
    grid(F01,'on'); set(F01,'XMinorGrid','on','YMinorGrid','on');  
    
    %----------% Chemical flux %----------%
    p(1,4).select(); F01=gca; 
    semilogx(F01,MFields.qc,zarray,'k','linewidth',linew); hold(F01,'on');
    xlabel(F01,'$\overline{Q}_c$','Fontsize',fontsize);
    set(F01,'ylim',[0 1],'Ytick',0:0.25:1)
    set(F01,'yaxislocation','right');       
    set(F01,'Box','on');
    grid(F01,'on'); set(F01,'XMinorGrid','on','YMinorGrid','on');    
 
    drawnow;
    
end

function plot_fluctuationsh(nfig,FFields,zarray,par,linew,fontsize)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT_FLUCATIONSH Plot the amplitude and phase of \hat{c} and \hat{phi}
%
%   n.b.: this function uses external-package panel-2.12
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    %----------%    Figure  Properties    %----------%    
    set(nfig,'position',[0 50 357 600]); 
    p=panel(nfig); p.pack(2,2); 
    p.fontsize=fontsize;
    p.title('Fluctuating state')
       
    p.margin    = [20 20 15 10]; % Margins outside [left bottom right top]
    p.de.margin = 11;            % Margins between sub-panels
    p(1).marginbottom = 20;    
    
    %----------% %----------%%----------%%----------%%----------%
    %----------%%----------% Amplitudes %----------%%----------%
    %----------%%----------%%----------%%----------%%----------%
    
    %----------% Amplitude of \hat{phi} %----------%
    p(1,1).select(); F01=gca; 
    plot(F01,log10(abs(FFields.phih)),zarray,'k','linewidth',linew);
    xlabel(F01,'$\vert\hat{\phi}\vert$','Fontsize',fontsize);
    ylabel(F01,'$z$'); 
    set(F01,'ylim',[0 1],'Ytick',0:0.25:1);
    set(F01,'yaxislocation','left');
    set(F01,'Box','on');
    grid(F01,'on'); set(F01,'XMinorGrid','on','YMinorGrid','on');     
    
    %----------% Amplitude of \hat{c} %----------%
    p(1,2).select(); F01=gca; 
    plot(F01,log10(abs(FFields.csh)),zarray,'k','linewidth',linew);
    xlabel(F01,'$\vert\hat{c}\vert$','Fontsize',fontsize);
    ylabel(F01,'$z$'); 
    set(F01,'ylim',[0 1],'Ytick',0:0.25:1);
    set(F01,'yaxislocation','right');
    set(F01,'Box','on');
    grid(F01,'on'); set(F01,'XMinorGrid','on','YMinorGrid','on');  
    

    %----------% %----------%%----------%%----------%%----------%
    %----------%%----------%    Phases  %----------%%----------%
    %----------%%----------%%----------%%----------%%----------%
    
 
    %----------% Phase-angle of \hat{phi} %----------%
    p(2,1).select(); F01=gca; 
    plot(F01,wrapTo2Pi(-angle(FFields.phih)),zarray,'k','linewidth',linew);
    xlabel(F01,'$\psi_{\hat{\phi}}$','Fontsize',fontsize);
    ylabel(F01,'$z$'); 
    set(F01,'xlim',[0 2*pi],'Xtick',[0 pi/2 pi 3*pi/2 2*pi],'Xticklabel',{'$0$','$\frac{\pi}{2}$','$\pi$','$\frac{3\pi}{2}$','$2\pi$'})
    set(F01,'ylim',[0 1],'Ytick',0:0.25:1);
    set(F01,'yaxislocation','left');
    set(F01,'Box','on');
    grid(F01,'on'); set(F01,'XMinorGrid','on','YMinorGrid','on');  
    
    %----------% Phase-angle of \hat{c} %----------%
    p(2,2).select(); F01=gca; 
    plot(F01,wrapTo2Pi(-angle(FFields.csh)),zarray,'k','linewidth',linew);
    xlabel(F01,'$\psi_{\hat{c_s}}$','Fontsize',fontsize);
    ylabel(F01,'$z$'); 
    set(F01,'xlim',[0 2*pi],'Xtick',[0 pi/2 pi 3*pi/2 2*pi],'Xticklabel',{'$0$','$\frac{\pi}{2}$','$\pi$','$\frac{3\pi}{2}$','$2\pi$'})
    set(F01,'ylim',[0 1],'Ytick',0:0.25:1);
    set(F01,'yaxislocation','right');
    set(F01,'Box','on');
    grid(F01,'on'); set(F01,'XMinorGrid','on','YMinorGrid','on');  
    
    
    drawnow;
    
end

function plot_fluctuations_ztspace(nfig,MFields,FFields,par,linew,fontsize,nperiods);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT_FLUCATIONS_ZTSPACE Plot the fluctuating variables \hat{phi},
%                         \hat{c}, \hat{f} and \hat{f}_c
%
%   n.b.: this function uses external-package panel-2.12
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %-----------%%----------% Defining arrays for z-t plots    %----------%%-----------%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    t = linspace(0,nperiods*par.tp/par.t0,par.ntime);  % Non-dimensional t-array
    z = linspace(0,1,par.nz);                          % Non-dimensional z-array
    [T,Z] = meshgrid(t,z); 

    phip  = par.delta0*real(exp(1i*par.omega*t')*FFields.phih);              % time-dependent fluctuating porosity
    clp   = par.delta0*real(exp(1i*par.omega*t')*FFields.csh/par.D_vol);     % time-dependent fluctuating liquid concentration
    fp    = par.delta0*real(exp(1i*par.omega*t')*FFields.qh);                % time-dependent fluctuating melt flux
    fcp   = par.delta0*real(exp(1i*par.omega*t')*FFields.qch);               % time-dependent fluctuating chemical flux
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %----------%%----------%        Options for figures    %----------%%----------%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    trange=linspace(0,nperiods*par.tp/par.t0,2*nperiods+1);
    phiprange2  = [-7 -4];
    clprange2   = [-2 2];
    fprange2    = [-4 -1];
    fcprange2   = [-3  0];
    
    nc=20;  % Number of sub-divisions in colormap
    cmap1 = colormap(summer(nc)); 
    cmap2 = colormap(autumn(nc)); cmap2 = cmap2(end:-1:1,:);
    cmap = [cmap1;cmap2]; 
    
    xtext2 = 0.40; ytext2=0.88; fonttext2 = 28; colortext=[0.0 0.0 0.0];
    x0_colbar1 = 0.10; y0_colbar1 = 0.05; xl_colbar1 = 0.35; yl_colbar1 = 0.03;
    x0_colbar2 = 0.57;  y0_colbar2 = 0.05; xl_colbar2 = 0.35; yl_colbar2 = 0.03; 

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %------------%%------------% Making figure %------------%%-----------%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    set(nfig,'position', [0, 100, 1000, 600]);
    
    p = panel(nfig);    
    p.pack({2/20 18/20}, 4);  
    p.fontsize=fontsize;
    
    %------------% Panel margins %------------%    

    p.de.margin=12;                                   %internal margins
    p(1,2).marginright = 20; p(2,2).marginright = 20; % Separate 2nd and 3rd columns
    p(1).marginbottom=8; p(2).margintop=0; 
    p.margin = [22 45 22 5];                          % margin[ = left bottom right top]

    
    %------------%%------------% FIRST ROW %------------%%------------%
    
    %------------% Plot SL and SL rate 1st column 
    p(1,1).select(); F01=gca;
    plot(F01,t,cos(par.omega*t),'--','linewidth',linew,'color',[0  0  0]); hold on; 
    plot(F01,t,sin(par.omega*t),':','linewidth',linew,'color',[0  0  0]); hold on; 
    plot(F01,t,(phip(:,end)-mean(phip(:,end)))./(max(phip(:,end))),'-','linewidth',linew,'color',[0 1 1]);
    set(F01,'Xlim',[trange(1) trange(end)],'Xtick',trange,'XTickLabel',[]);
    set(F01,'Ylim',[-1 1],'Ytick',[-1 0 1],'YTickLabel',{'$-1$' '0' '$1$'});
    set(F01,'Box','on');
    hold(F01,'on');
    grid(F01,'on');  

   
    %------------% Plot SL and SL rate 2nd column     
    p(1,2).select(); F01=gca;
    
    plot(F01,t,cos(par.omega*t),'--','linewidth',linew,'color',[0  0  0]); hold on; 
    plot(F01,t,sin(par.omega*t),':','linewidth',linew,'color',[0  0  0]); hold on; 
    plot(F01,t,(clp(:,end)-mean(clp(:,end)))./(max(clp(:,end))),'-','linewidth',linew,'color',[1 0 0]);
    set(F01,'Xlim',[trange(1) trange(end)],'Xtick',trange,'XTickLabel',[]);
    set(F01,'Ylim',[-1 1],'Ytick',[-1 0 1],'YTickLabel',[]);
    if (strcmp(par.Gammap,'off')==1)
        set(F01,'YTickLabel',[]);
    else
        set(F01,'YTickLabel',[],'YAxisLocation','right');
    end
    set(F01,'Box','on');
    hold(F01,'on');
    grid(F01,'on');  

    %------------% Plot SL and SL rate 3rd column     
    p(1,3).select(); F01=gca;
    
    plot(F01,t,cos(par.omega*t),'--','linewidth',linew,'color',[0  0  0]); hold on; 
    plot(F01,t,sin(par.omega*t),':','linewidth',linew,'color',[0  0  0]); hold on; 
    plot(F01,t,(fp(:,end)-mean(fp(:,end)))./(max(fp(:,end))),'-','linewidth',linew,'color',[0 0 1]);
    set(F01,'Xlim',[trange(1) trange(end)],'Xtick',trange,'XTickLabel',[]);
    set(F01,'Ylim',[-1 1],'Ytick',[-1 0 1],'YTickLabel',[]);
    if (strcmp(par.Gammap,'off')==1)
        set(F01,'YTickLabel',[]);
    else
        set(F01,'YTickLabel',[],'YAxisLocation','right');
    end
    set(F01,'Box','on');
    hold(F01,'on');
    grid(F01,'on');  

    
    %------------% Plot SL and SL rate 4th column     
    p(1,4).select(); F01=gca;
    
    plot(F01,t,cos(par.omega*t),'--','linewidth',linew,'color',[0  0  0]); hold on; 
    plot(F01,t,sin(par.omega*t),':','linewidth',linew,'color',[0  0  0]); hold on; 
    plot(F01,t,(fcp(:,end)-mean(fcp(:,end)))./(max(fcp(:,end))),'-','linewidth',linew,'color',[0.5 0 0.5]);
    set(F01,'Xlim',[trange(1) trange(end)],'Xtick',trange,'XTickLabel',[]);
    set(F01,'Ylim',[-1 1],'Ytick',[-1 0 1],'YTickLabel',{'$-1$' '0' '$1$'},'YAxisLocation','right');
    set(F01,'Box','on');
    hold(F01,'on');
    grid(F01,'on'); 
    

    %------------%%------------%  PLOT POROSITY  %------------%%------------%
    
    [phip_m,logscale,ticks_rescaled] = log_negative(phip,phiprange2); % Rescale field to display log of negative values
    str_colorbar_label = sprintf('$\\phi^{\\prime}$');                % Label of field 
    
    p(2,1).select(); F01=gca;   
        
    [FC,h]=contourf(T,Z,phip_m(1:end,:)',linspace(ticks_rescaled(1),ticks_rescaled(end),2*nc+1)); h.LineStyle='none'; hold on;
    caxis(F01,[min(ticks_rescaled) max(ticks_rescaled)]); 
    colormap(F01,cmap);
         
    %------------% Managing Axis
    set(F01,'Box','on');
    set(F01,'Xlim',[trange(1) trange(end)],'Xtick',trange,'XTickLabel',{'0' '$\frac{t_p}{2}$' '$t_p$' '$\frac{3 t_p}{2}$' '$2t_p$'});
    F01.XGrid = 'on'; F01.YGrid = 'on'; F01.Layer='top';
    text(-0.3,0.5,'$z$','Units','normalized','fontsize',24);
    text(xtext2,ytext2,str_colorbar_label,'Units','normalized','fontsize',fonttext2,'color',colortext); 
        
    %------------% Managing colorbar
    colbar = draw_colorbar(F01,caxis,logscale,[x0_colbar1 y0_colbar1 xl_colbar1 yl_colbar1],[-0.12 1.5],str_colorbar_label);
    colbar.YAxisLocation = 'top';  
    
    
    %------------%%------------%  PLOT LIQUID CONCENTRATION  %------------%%------------%
    
    [clp_m,logscale,ticks_rescaled] = log_negative(clp,clprange2); % Rescale field to display log of negative values
    str_colorbar_label = sprintf('$c_l^{\\prime}$');              % Label of field 
    
     p(2,2).select(); F01=gca;
    
    [FC,h]=contourf(T,Z,clp_m(1:end,:)',linspace(ticks_rescaled(1),ticks_rescaled(end),2*nc+1)); h.LineStyle='none'; hold on;
    caxis(F01,[min(ticks_rescaled) max(ticks_rescaled)]); 
    colormap(F01,cmap);
    
    %------------% Managing Axis
    set(F01,'Box','on');
    set(F01,'Xlim',[trange(1) trange(end)],'Xtick',trange,'XTickLabel',{'0' '$\frac{t_p}{2}$' '$t_p$' '$\frac{3 t_p}{2}$' '$2t_p$'});
    set(F01,'YTickLabel',[]);
    set(F01,'YAxisLocation','right');
    F01.XGrid = 'on'; F01.YGrid = 'on'; F01.Layer='top';
    text(xtext2,ytext2,str_colorbar_label,'Units','normalized','fontsize',fonttext2,'color',colortext); 

    %------------% Managing colorbar
    colbar = draw_colorbar(F01,caxis,logscale,[x0_colbar1 y0_colbar1 xl_colbar1 yl_colbar1],[-0.12 -0.2],str_colorbar_label);
    
    
    %------------%%------------%  PLOT MELT FLUX  %------------%%------------%
    
    [fp_m,logscale,ticks_rescaled] = log_negative(fp,fprange2); % Rescale field to display log of negative values
    str_colorbar_label = sprintf('$Q^{\\prime}$');              % Label of field  
    
    p(2,3).select(); F01=gca;

    [FC,h]=contourf(T,Z,fp_m(1:end,:)',linspace(ticks_rescaled(1),ticks_rescaled(end),2*nc+1)); h.LineStyle='none'; hold on;
    caxis(F01,[min(ticks_rescaled) max(ticks_rescaled)]); 
    colormap(F01,cmap);

    %------------% Managing Axis
    set(F01,'Box','on');
    set(F01,'Xlim',[trange(1) trange(end)],'Xtick',trange,'XTickLabel',{'0' '$\frac{t_p}{2}$' '$t_p$' '$\frac{3 t_p}{2}$' '$2t_p$'});
    set(F01,'YTickLabel',[]);
    F01.XGrid = 'on'; F01.YGrid = 'on'; F01.Layer='top';
    text(xtext2,ytext2,str_colorbar_label,'Units','normalized','fontsize',fonttext2,'color',colortext); 

    %------------% Managing colorbar
    colbar = draw_colorbar(F01,caxis,logscale,[x0_colbar2 y0_colbar2 xl_colbar2 yl_colbar2],[-0.12 1.5],str_colorbar_label);
    colbar.YAxisLocation = 'top';     

    
    %------------%%------------%  PLOT CARBON FLUX  %------------%%------------%
    
    [fcp_m,logscale,ticks_rescaled] = log_negative(fcp,fcprange2); % Rescale field to display log of negative values
    str_colorbar_label = sprintf('$Q_c^{\\prime}$');               % Label of field  
   
    p(2,4).select(); F01=gca;

    [FC,h]=contourf(T,Z,fcp_m(1:end,:)',linspace(ticks_rescaled(1),ticks_rescaled(end),2*nc+1)); h.LineStyle='none'; hold on;
    caxis(F01,[min(ticks_rescaled) max(ticks_rescaled)]); 
    colormap(F01,cmap);

    %------------% Managing Axis
    set(F01,'Box','on');
    text(1.2,0.5,'$z$','Units','normalized','fontsize',24);
    set(F01,'Xlim',[trange(1) trange(end)],'Xtick',trange,'XTickLabel',{'0' '$\frac{t_p}{2}$' '$t_p$' '$\frac{3 t_p}{2}$' '$2t_p$'});
    set(F01,'YAxisLocation','right');
    F01.XGrid = 'on'; F01.YGrid = 'on'; F01.Layer='top';
    text(xtext2,ytext2,str_colorbar_label,'Units','normalized','fontsize',fonttext2,'color',colortext); 

    %------------% Managing colorbar
    colbar = draw_colorbar(F01,caxis,logscale,[x0_colbar2 y0_colbar2 xl_colbar2 yl_colbar2],[-0.12 -0.2],str_colorbar_label);
        
    
    
    drawnow;
    
    
        function [field_m,logscale,ticks_rescaled] = log_negative(field,range)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %  LOG_NEGATIVE rescales fields to plot log of negative values
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            field0=field;
            field(abs(field)<10^range(1))=10^range(1)*sign(field(abs(field)<10^range(1)));
            field(abs(field)>10^range(2))=10^range(2)*sign(field(abs(field)>10^range(2)));

            % Defining bounds
            logscale.min_sat  = range(1);
            logscale.max_sat  = range(2); 
            logscale.min      = min(min(log10(abs(field0))));
            logscale.max      = max(max(log10(abs(field0))));
            logscale.diff_sat = logscale.max_sat - logscale.min_sat + 1; % +1 is to take into account lower bound for saturation
            logscale.end_sat  = round(logscale.min_sat - logscale.diff_sat);

            % Rescaling array
            field_m = zeros(size(field));
            field_m(field>0)  = log10(field(field>0));
            field_m(field<0)  = logscale.end_sat  - log10(-field(field<0)) + logscale.max_sat;
            field_m(field==0) = logscale.min_sat; 
            field_m(field_m<logscale.end_sat)=logscale.end_sat;

            ticks_rescaled =  logscale.end_sat:1:range(2);    


        end
    
        function colbar = draw_colorbar(fig,caxis,logscale,pos_colorbar,pos_label,str_colorbar_label)
           
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %  DRAW_COLORBAR displays colorbar and manages label/ticks-label appearance
        %    Output 
        %       colbar : colorbar properties object
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            colbar=colorbar(fig,'Location','southoutside','Position',pos_colorbar,'Units','normalized'); % colorbar position
            colbar.Limits= caxis;
            colbar.Ticks = ticks_rescaled;
            colbar.TickDirection='out';  colbar.TickLength=0.02;
            colbar.TickLabelInterpreter='LaTex'; 
            for I=1:numel(colbar.TickLabels)
                n=str2num(colbar.TickLabels{I});
                if (n<logscale.min_sat)
                    colbar.TickLabels{I} = ['$-10^{',num2str(logscale.max_sat+(logscale.end_sat-n)),'}$'];
                else
                    colbar.TickLabels{I} = ['$10^{',num2str(n),'}$'];
                end
            end  
            colbar.Label.Interpreter='LaTex'; colbar.Label.String=str_colorbar_label; % colorbar label
            colbar.Label.Units='normalized';  colbar.Label.Position=pos_label;        % position of label relative to colorbar
            colbar.FontSize=15; colbar.Label.FontSize=18;

        end
end

