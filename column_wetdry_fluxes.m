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
%----------% Reading input parameters %----------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

par = input_parameters();

zarray = linspace(0,1,par.nz);

%----------% Plotting parameters %----------%
linew = 3;  % linewidth
fonts = 18; % fontsize

%------------% Range of fluctuating variables in plots %------------%
frange        = [-4 -1]; % melt flux
fcrange       = [-3  0]; % carbon flux

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----------%    Compute solution    %----------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%----------%  Defining model parameters array %----------%
Harray   = [par.H  par.Hdry   par.Hdry];
RHSarray = {'on'  'on'    'off'};
epsM     = 1e-8; % to avoid singularity for M = 0

tparray = [23 41 100]; % tp-array wjTh tp in kyr

%----------%  Loop over model parameters  %----------%
par.verb = "off";
for iH = 1:length(Harray)   
    
    %----------% Updating parameters %----------%
    par.G = par.Fmax*Harray(iH)/par.Hdry;
    par.M = par.Fmax*(Harray(iH)/par.Hdry - 1) + epsM; 
    par.Gammap = RHSarray{iH};
    
    %----------% Calculate mean variables %----------%
    [~,~,MFields.cs,MFields.phi,~,~] = mean_analytical(zarray,par); %calculate base state analytically
    
    for jT = 1:length(tparray);
        
        fprintf('\n Model %3d.%1d : H = %4.1f km, Tp = %4.1f kyr',iH,jT,Harray(iH),tparray(jT))
        
        par.tp    = tparray(jT)*1e3; 
        par.omega = 2.*pi/par.tp*par.t0;
        
        %----------% Get other mean variables %----------%
        [MFields.W,MFields.q,MFields.qc] = get_other_mfields(MFields.phi,MFields.cs,par);

        %----------% Calculate fluctuating variables %----------%
        [~,~,FFields.csh,FFields.phiH] = fluctuations(zarray,par);
        %----------% Get other fluctuating variables %----------%
        [FFields.Wh,FFields.qh,FFields.qch] = get_other_ffields(MFields.phi,MFields.cs,FFields.phiH,FFields.csh,par); 

        array_MFields{iH,jT} = MFields; 
        array_FFields{iH,jT} = FFields;
        array_par{iH,jT} = par; 
    
    end
    
end
fprintf('\n... end of loop over models \n\n')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----------%    Making figures    %----------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nfig=0;

%----------% Melt and carbon flux in wet model %----------%
nfig = nfig+1; figure(nfig); 
plot_wet_fluxes(nfig,array_MFields,array_FFields,array_par,fonts,linew,frange,fcrange);

%----------% Melt flux in dry model %----------%
nfig = nfig+1; figure(nfig); 
plot_dry_meltflux(nfig,array_MFields,array_FFields,array_par,fonts,linew,frange);

%----------% Carbon flux in basal-flux model %----------%
nfig = nfig+1; figure(nfig); 
plot_basal_carbonflux(nfig,array_MFields,array_FFields,array_par,fonts,linew,fcrange);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%         LOCAL FUNCTIONS         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_wet_fluxes(nfig,array_MFields,array_FFields,array_par,fontsize,linew,frange,fcrange); 


    %------------% Figure properties %------------%
    set(nfig,'position', [0, 100, 700, 1000]);

    % Create panel frame
    p = panel(nfig);

    % Create panel 2x1
    p.pack({10/210 100/210 100/210 }, 3);  %
    
    p.fontsize=fontsize;

    %------------% Other panel properties %------------%    

    p.de.margin=8; %internal margins
    p(1).marginbottom=5; p(2).margintop=0; 
    p(2).marginbottom=5; p(3).margintop=0; % Tighten 2nd and 3rd rows

    p.margin = [20 40 15 9];   %% margin[left bottom right top]

    
    nc=20; % number of colors in colorbar
    cmap1 = colormap(summer(nc)); 
    cmap2 = colormap(autumn(nc)); cmap2 = cmap2(end:-1:1,:);
    cmap = [cmap1;cmap2]; 
    
    zrange2_lin = [0.0 1]; 
    zticks2_lin = [0 0.25 0.5 0.75 1.0]; zticks2_lin_label = {'0' '0.25','0.5','0.75','1.0'};

    alphagrid=0.2;
    nperiods = 2.0; 

    slcolor = [0.5  0.5  0.5];
    lstyle  = {'-', '--', ':'};

    x0_colbar1 = 0.23; y0_colbar1 = 0.03; xl_colbar1 = 0.6; yl_colbar1 = 0.015;
    
    
    iH = 1;
    
    
    %------------%%------------%  PLOT MELT FLUX  %------------%%------------%
    str_colorbar_label = sprintf('$Q^{\\prime}$');
    fprintf('\n PLOTTING melt flux in wet models... ')
    
    %------------% Loop over models %------------%
    nTps    = size(array_MFields,2);
    for jT = 1:nTps
        
        %------------% Load arrays %------------%
        par      = array_par{1,jT}; 
        FFields  = array_FFields{1,jT}; 

        t = linspace(0,nperiods*par.tp/par.t0,par.ntime);
        trange=linspace(0,nperiods,5);    
        t_rescale = par.t0/par.tp;
        z = linspace(0,1,par.nz); 

        fp          = par.delta0*real(exp(1i*par.omega*t')*FFields.qh);         

        %zprime_dry = 1 - par.Hdry/par.H;
        [T,Z] = meshgrid(t*t_rescale,z); 

        %------------% Plot SL and SL rate %------------%
        p(1,jT).select(); F_row1=gca;
        if (iH == 1)
            p(iH,jT).select(); 
            plot(F_row1,t*t_rescale,cos(par.omega*t),'--','linewidth',linew,'color',slcolor); hold on;
            plot(F_row1,t*t_rescale,sin(par.omega*t),':','linewidth',linew,'color',slcolor); hold on; 
            set(F_row1,'Ylim',[-1 1],'Ytick',[-1 1]);
            set(F_row1,'YTickLabel',{'$-1$' '$+1$'});
            title(F_row1,sprintf('$t_p$=%3d kyr',par.tp*1e-3)); 
            if (jT==2)
                set(F_row1,'YTickLabel',[]);
            elseif (jT==3) 
                set(F_row1,'YAxisLocation','right');
            end
            set(F_row1,'Box','on');
            grid(F_row1,'on');  F_row1.GridAlpha = alphagrid;  
        end
        plot(F_row1,t*t_rescale,(fp(:,end)-mean(fp(:,end)))./(max(fp(:,end))),'linewidth',linew,'linestyle',lstyle{iH},'color',[0.0 0.0 1.0]); hold on;
        set(F_row1,'Xlim',[0 2],'Xtick',[0 0.5 1.0 3/2 2],'XTickLabel',[]);

        %------------%%------------% Plotting flux %------------%%------------%
        p(1+iH,jT).select(); F01=gca;

        [fp_m,logscale,ticks_rescaled] = log_negative(fp,frange); % Rescale field to display log of negative values    

        [FC,h]=contourf(F01,T,Z,fp_m(1:end,:)',linspace(ticks_rescaled(1),ticks_rescaled(end),2*nc+1)); h.LineStyle='none'; hold on;
        caxis(F01,[min(ticks_rescaled) max(ticks_rescaled)]); 
        colormap(F01,cmap);

        %------------% Managing Axis
        set(F01,'Box','on');
        set(F01,'Xlim',[trange(1) trange(end)],'Xtick',trange);
        %plot(F01,[t(1)*t_rescale t(end)*t_rescale],[zprime_dry zprime_dry],'k--','linewidth',1.5);
        set(F01,'XTickLabel',[]);
        set(F01,'Ylim',zrange2_lin);
        set(F01,'YTick',zticks2_lin);
        if (jT == 1)
           set(F01,'YTickLabel',zticks2_lin_label);
           text(-0.25,0.5,'$z$','Units','normalized','fontsize',fontsize);
        elseif (jT == 2) 
           set(F01,'YTickLabel',[]);
        elseif (jT ==3)
           set(F01,'YTickLabel',zticks2_lin_label);
           set(F01,'YAxisLocation','right');
        end
        F01.XGrid = 'on'; F01.YGrid = 'on'; F01.GridAlpha = alphagrid; F01.Layer='top';
        if (jT == 2)          
            %------------% Managing ColorBar
            colbar = draw_colorbar(F01,caxis,logscale,[x0_colbar1 y0_colbar1 xl_colbar1 yl_colbar1],[-0.1 1.5],str_colorbar_label,ticks_rescaled);
            colbar.YAxisLocation = 'top'; 
        else
        end   
    
    end
    drawnow;
    fprintf(' ... DONE \n\n');
  
     
   
    %------------%%------------%  PLOT CARBON FLUX  %------------%%------------%
    str_colorbar_label = sprintf('$Q_c^{\\prime}$');
    fprintf('\n PLOTTING carbon flux in wet models... ')
    
     %------------% Loop over models %------------%
    nTps    = size(array_MFields,2);
    for jT = 1:nTps
        
        %------------% Load arrays %------------%
        par      = array_par{1,jT}; 
        FFields  = array_FFields{1,jT}; 

        t = linspace(0,nperiods*par.tp/par.t0,par.ntime);
        trange=linspace(0,nperiods,5);    
        t_rescale = par.t0/par.tp;
        z = linspace(0,1,par.nz); 

        fp          = par.delta0*real(exp(1i*par.omega*t')*FFields.qch); 

        zprime_dry = 1 - par.Hdry/par.H;
        [T,Z] = meshgrid(t*t_rescale,z); 

        %------------% Plot SL and SL rate %------------%
        p(1,jT).select(); F_row1=gca;
        if (iH == 1)
            p(iH,jT).select(); 
            plot(F_row1,t*t_rescale,cos(par.omega*t),'--','linewidth',linew,'color',slcolor); hold on;
            plot(F_row1,t*t_rescale,sin(par.omega*t),':','linewidth',linew,'color',slcolor); hold on; 
            set(F_row1,'Ylim',[-1 1],'Ytick',[-1 1]);
            set(F_row1,'YTickLabel',{'$-1$' '$+1$'});
            title(F_row1,sprintf('$t_p$=%3d kyr',par.tp*1e-3)); 
            if (jT==2)
                set(F_row1,'YTickLabel',[]);
            elseif (jT==3) 
                set(F_row1,'YAxisLocation','right');
            end
            set(F_row1,'Box','on');
            grid(F_row1,'on');  F_row1.GridAlpha = alphagrid;  
        end
        plot(F_row1,t*t_rescale,(fp(:,end)-mean(fp(:,end)))./(max(fp(:,end))),'linewidth',linew,'linestyle',lstyle{iH},'color',[0.5 0.0 0.5]); hold on;
        set(F_row1,'Xlim',[0 2],'Xtick',[0 0.5 1.0 3/2 2],'XTickLabel',[]);

        %------------%%------------% Plotting carbon flux %------------%%------------%
        p(2+iH,jT).select(); F01=gca;

        [fp_m,logscale,ticks_rescaled] = log_negative(fp,fcrange); % Rescale field to display log of negative values    

        [FC,h]=contourf(F01,T,Z,fp_m(1:end,:)',linspace(ticks_rescaled(1),ticks_rescaled(end),2*nc+1)); h.LineStyle='none'; hold on;
        caxis(F01,[min(ticks_rescaled) max(ticks_rescaled)]); 
        colormap(F01,cmap);

        %------------% Managing Axis
        set(F01,'Box','on');
        set(F01,'Xlim',[trange(1) trange(end)],'Xtick',trange);
        %plot(F01,[t(1)*t_rescale t(end)*t_rescale],[zprime_dry zprime_dry],'k--','linewidth',1.5);
        set(F01,'XTickLabel',{'0' '$\frac{t_p}{2}$' '$t_p$' '$\frac{3 t_p}{2}$' '$2t_p$'});
        set(F01,'Ylim',zrange2_lin);
        set(F01,'YTick',zticks2_lin);
        if (jT == 1)
           set(F01,'YTickLabel',zticks2_lin_label);
           text(-0.25,0.5,'$z$','Units','normalized','fontsize',fontsize);
        elseif (jT == 2) 
           set(F01,'YTickLabel',[]);
        elseif (jT ==3)
           set(F01,'YTickLabel',zticks2_lin_label);
           set(F01,'YAxisLocation','right');
        end
        F01.XGrid = 'on'; F01.YGrid = 'on'; F01.GridAlpha = alphagrid; F01.Layer='top';
        if (jT == 2)          
            %------------% Managing ColorBar
            colbar = draw_colorbar(F01,caxis,logscale,[x0_colbar1 y0_colbar1 xl_colbar1 yl_colbar1],[-0.1 -0.2],str_colorbar_label,ticks_rescaled);
        else
        end   
    
    end
    drawnow;
    fprintf(' ... DONE \n\n');

    
end

function plot_dry_meltflux(nfig,array_MFields,array_FFields,array_par,fontsize,linew,frange); 


    %------------% Figure properties %------------%
    set(nfig,'position', [0, 100, 700, 400]);

    % Create panel frame
    p = panel(nfig);

    % Create panel 2x1
    p.pack({10/110 100/110}, 3);  %
    
    p.fontsize=fontsize;

    %------------% Other panel properties %------------%    

    p.de.margin=8; %internal margins
    p(1).marginbottom=5; p(2).margintop=0; 

    p.margin = [20 30 15 9];   %% margin[left bottom right top]

    nc=20; % number of colors in colorbar
    cmap1 = colormap(summer(nc)); 
    cmap2 = colormap(autumn(nc)); cmap2 = cmap2(end:-1:1,:);
    cmap = [cmap1;cmap2]; 
    
    zrange2_lin = [0.0 1]; 
    zticks2_lin = [0 0.25 0.5 0.75 1.0]; zticks2_lin_label = {'0' '0.25','0.5','0.75','1.0'};

    alphagrid=0.2;
    nperiods = 2.0; 

    slcolor = [0.5  0.5  0.5];
    lstyle  = {'-', '--', ':'};

    x0_colbar1 = 0.23; y0_colbar1 = 0.08; xl_colbar1 = 0.6; yl_colbar1 = 0.02;
    
    iH = 2;
    
    %------------%%------------%  PLOT MELT FLUX  %------------%%------------%
    str_colorbar_label = sprintf('$Q^{\\prime}$');
    fprintf('\n PLOTTING melt flux in dry models... ')
    
    %------------% Loop over models %------------%
    nTps    = size(array_MFields,2);
    for jT = 1:nTps
        

        %------------% Load arrays %------------%
        par      = array_par{iH,jT}; 
        FFields  = array_FFields{iH,jT}; 

        t = linspace(0,nperiods*par.tp/par.t0,par.ntime);
        trange=linspace(0,nperiods,5);    
        t_rescale = par.t0/par.tp;
        z = linspace(0,1,par.nz); 

        fp          = par.delta0*real(exp(1i*par.omega*t')*FFields.qh);         

        [T,Z] = meshgrid(t*t_rescale,z); 

        %------------% Plot SL and SL rate %------------%
        p(1,jT).select(); F_row1=gca;

            plot(F_row1,t*t_rescale,cos(par.omega*t),'--','linewidth',linew,'color',slcolor); hold on;
            plot(F_row1,t*t_rescale,sin(par.omega*t),':','linewidth',linew,'color',slcolor); hold on; 
            set(F_row1,'Ylim',[-1 1],'Ytick',[-1 1]);
            set(F_row1,'YTickLabel',{'$-1$' '$+1$'});
            title(F_row1,sprintf('$t_p$=%3d kyr',par.tp*1e-3)); 
            if (jT==2)
                set(F_row1,'YTickLabel',[]);
            elseif (jT==3) 
                set(F_row1,'YAxisLocation','right');
            end
            set(F_row1,'Box','on');
            grid(F_row1,'on');  F_row1.GridAlpha = alphagrid;  

        plot(F_row1,t*t_rescale,(fp(:,end)-mean(fp(:,end)))./(max(fp(:,end))),'linewidth',linew,'linestyle',lstyle{iH},'color',[0.0 0.0 1.0]); hold on;
        set(F_row1,'Xlim',[0 2],'Xtick',[0 0.5 1.0 3/2 2],'XTickLabel',[]);

        %------------%%------------% Plotting flux %------------%%------------%
        p(2,jT).select(); F01=gca;

        [fp_m,logscale,ticks_rescaled] = log_negative(fp,frange); % Rescale field to display log of negative values    

        [FC,h]=contourf(F01,T,Z,fp_m(1:end,:)',linspace(ticks_rescaled(1),ticks_rescaled(end),2*nc+1)); h.LineStyle='none'; hold on;
        caxis(F01,[min(ticks_rescaled) max(ticks_rescaled)]); 
        colormap(F01,cmap);

        %------------% Managing Axis
        set(F01,'Box','on');
        set(F01,'Xlim',[trange(1) trange(end)],'Xtick',trange);
        set(F01,'XTickLabel',{'0' '$\frac{t_p}{2}$' '$t_p$' '$\frac{3 t_p}{2}$' '$2t_p$'});
        set(F01,'Ylim',zrange2_lin);
        set(F01,'YTick',zticks2_lin);
        if (jT == 1)
           set(F01,'YTickLabel',zticks2_lin_label);
           text(-0.25,0.5,'$z$','Units','normalized','fontsize',fontsize);
        elseif (jT == 2) 
           set(F01,'YTickLabel',[]);
        elseif (jT ==3)
           set(F01,'YTickLabel',zticks2_lin_label);
           set(F01,'YAxisLocation','right');
        end
        F01.XGrid = 'on'; F01.YGrid = 'on'; F01.GridAlpha = alphagrid; F01.Layer='top';
        if (jT == 2)          
            %------------% Managing ColorBar
            colbar = draw_colorbar(F01,caxis,logscale,[x0_colbar1 y0_colbar1 xl_colbar1 yl_colbar1],[-0.1 -.2],str_colorbar_label,ticks_rescaled);
        else
        end   
    
    end
    drawnow;
    fprintf(' ... DONE \n\n');
  
    
end

function plot_basal_carbonflux(nfig,array_MFields,array_FFields,array_par,fontsize,linew,frange); 


    %------------% Figure properties %------------%
    set(nfig,'position', [0, 100, 700, 400]);

    % Create panel frame
    p = panel(nfig);

    % Create panel 2x1
    p.pack({10/110 100/110}, 3);  %
    
    p.fontsize=fontsize;

    %------------% Other panel properties %------------%    

    p.de.margin=8; %internal margins
    p(1).marginbottom=5; p(2).margintop=0; 

    p.margin = [20 30 15 9];   %% margin[left bottom right top]

    nc=20; % number of colors in colorbar
    cmap1 = colormap(summer(nc)); 
    cmap2 = colormap(autumn(nc)); cmap2 = cmap2(end:-1:1,:);
    cmap = [cmap1;cmap2]; 
    
    zrange2_lin = [0.0 1]; 
    zticks2_lin = [0 0.25 0.5 0.75 1.0]; zticks2_lin_label = {'0' '0.25','0.5','0.75','1.0'};

    alphagrid=0.2;
    nperiods = 2.0; 

    slcolor = [0.5  0.5  0.5];
    lstyle  = {'-', '--', ':'};

    x0_colbar1 = 0.23; y0_colbar1 = 0.08; xl_colbar1 = 0.6; yl_colbar1 = 0.02;
    
    iH = 3;
    
    %------------%%------------%  PLOT CARBON FLUX  %------------%%------------%
    str_colorbar_label = sprintf('$Q_c^{\\prime}$');
    fprintf('\n PLOTTING carbon flux in basal-flux models... ')
    
    %------------% Loop over models %------------%
    nTps    = size(array_MFields,2);
    for jT = 1:nTps
        

        %------------% Load arrays %------------%
        par      = array_par{iH,jT}; 
        FFields  = array_FFields{iH,jT}; 

        t = linspace(0,nperiods*par.tp/par.t0,par.ntime);
        trange=linspace(0,nperiods,5);    
        t_rescale = par.t0/par.tp;
        z = linspace(0,1,par.nz); 

        fp    = par.delta0*real(exp(1i*par.omega*t')*FFields.qch);         

        [T,Z] = meshgrid(t*t_rescale,z); 

        %------------% Plot SL and SL rate %------------%
        p(1,jT).select(); F_row1=gca;

            plot(F_row1,t*t_rescale,cos(par.omega*t),'--','linewidth',linew,'color',slcolor); hold on;
            plot(F_row1,t*t_rescale,sin(par.omega*t),':','linewidth',linew,'color',slcolor); hold on; 
            set(F_row1,'Ylim',[-1 1],'Ytick',[-1 1]);
            set(F_row1,'YTickLabel',{'$-1$' '$+1$'});
            title(F_row1,sprintf('$t_p$=%3d kyr',par.tp*1e-3)); 
            if (jT==2)
                set(F_row1,'YTickLabel',[]);
            elseif (jT==3) 
                set(F_row1,'YAxisLocation','right');
            end
            set(F_row1,'Box','on');
            grid(F_row1,'on');  F_row1.GridAlpha = alphagrid;  

        plot(F_row1,t*t_rescale,(fp(:,end)-mean(fp(:,end)))./(max(fp(:,end))),'linewidth',linew,'linestyle',lstyle{iH},'color',[0.5 0.0 0.5]); hold on;
        set(F_row1,'Xlim',[0 2],'Xtick',[0 0.5 1.0 3/2 2],'XTickLabel',[]);

        %------------%%------------% Plotting flux %------------%%------------%
        p(2,jT).select(); F01=gca;

        [fp_m,logscale,ticks_rescaled] = log_negative(fp,frange); % Rescale field to display log of negative values    

        [FC,h]=contourf(F01,T,Z,fp_m(1:end,:)',linspace(ticks_rescaled(1),ticks_rescaled(end),2*nc+1)); h.LineStyle='none'; hold on;
        caxis(F01,[min(ticks_rescaled) max(ticks_rescaled)]); 
        colormap(F01,cmap);

        %------------% Managing Axis
        set(F01,'Box','on');
        set(F01,'Xlim',[trange(1) trange(end)],'Xtick',trange);
        %plot(F01,[t(1)*t_rescale t(end)*t_rescale],[zprime_dry zprime_dry],'k--','linewidth',1.5);
        set(F01,'XTickLabel',{'0' '$\frac{t_p}{2}$' '$t_p$' '$\frac{3 t_p}{2}$' '$2t_p$'});
        set(F01,'Ylim',zrange2_lin);
        set(F01,'YTick',zticks2_lin);
        if (jT == 1)
           set(F01,'YTickLabel',zticks2_lin_label);
           text(-0.25,0.5,'$z$','Units','normalized','fontsize',fontsize);
        elseif (jT == 2) 
           set(F01,'YTickLabel',[]);
        elseif (jT ==3)
           set(F01,'YTickLabel',zticks2_lin_label);
           set(F01,'YAxisLocation','right');
        end
        F01.XGrid = 'on'; F01.YGrid = 'on'; F01.GridAlpha = alphagrid; F01.Layer='top';
        if (jT == 2)          
            %------------% Managing ColorBar
            colbar = draw_colorbar(F01,caxis,logscale,[x0_colbar1 y0_colbar1 xl_colbar1 yl_colbar1],[-0.1 -.2],str_colorbar_label,ticks_rescaled);
        else
        end   
    
    end
    drawnow;
    fprintf(' ... DONE \n\n');
  
    
end

function [field_m,logscale,ticks_rescaled] = log_negative(field,range)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  LOG_NEGATIVE rescales fields to plot log of negative values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    field0=field;
    field(abs(field)<10^range(1))=10^range(1)*sign(field(abs(field)<10^range(1)));
    field(abs(field)>10^range(2))=10^range(2)*sign(field(abs(field)>10^range(2)));

    %------------% Defining bounds
    logscale.min_sat  = range(1);
    logscale.max_sat  = range(2); 
    logscale.min      = min(min(log10(abs(field0))));
    logscale.max      = max(max(log10(abs(field0))));
    logscale.diff_sat = logscale.max_sat - logscale.min_sat + 1; %% +1 is to take into account lower bound for saturation
    logscale.end_sat  = round(logscale.min_sat - logscale.diff_sat);

    %------------% Rescaling array
    field_m = zeros(size(field));
    field_m(field>0)  = log10(field(field>0));
    field_m(field<0)  = logscale.end_sat  - log10(-field(field<0)) + logscale.max_sat;
    field_m(field==0) = logscale.min_sat; 
    field_m(field_m<logscale.end_sat)=logscale.end_sat;

    ticks_rescaled =  logscale.end_sat:1:range(2);    


end

function colbar = draw_colorbar(fig,caxis,logscale,pos_colorbar,pos_label,str_colorbar_label,ticks_rescaled)

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

