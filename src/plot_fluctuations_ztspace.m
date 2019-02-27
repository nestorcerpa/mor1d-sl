function plot_fluctuations_ztspace(nfig,MFields,FFields,par,linew,fontsize,nperiods);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %------------%%------------% Defining arrays for z-t plots %------------%%------------%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    t = linspace(0,nperiods*par.tp/par.t0,par.ntime);  % Non-dimensional t-array
    z = linspace(0,1,par.nz);                          % Non-dimensional z-array
    [T,Z] = meshgrid(t,z); 

    phip  = par.delta0*real(exp(1i*par.omega*t')*FFields.phih);              % time-dependent fluctuating porosity
    clp   = par.delta0*real(exp(1i*par.omega*t')*FFields.csh/par.D_vol);     % time-dependent fluctuating liquid concentration
    fp    = par.delta0*real(exp(1i*par.omega*t')*FFields.fh);                % time-dependent fluctuating melt flux
    fcp   = par.delta0*real(exp(1i*par.omega*t')*FFields.fch);               % time-dependent fluctuating chemical flux

    slev   = real(exp(1i*par.omega*t));                           % sea-level
    slrate = real(1i*par.omega*exp(1i*par.omega*t));              % rate of sea-level
    
    
    %----------% Options for figures %----------%
    trange=linspace(0,nperiods*par.tp/par.t0,2*nperiods+1);
    phiprange2  = [-7 -4];
    clprange2   = [-2 2];
    fprange2    = [-4 -1];
    fcprange2   = [-3  0];
    
    nc=20;  % Number of sub-divisions in colormap
    cmap1 = colormap(summer(nc)); %cmap1 = cmap1(end:-1:1,:);
    cmap2 = colormap(autumn(nc)); cmap2 = cmap2(end:-1:1,:);
    cmap = [cmap1;cmap2]; 
    
    xtext2 = 0.40; ytext2=0.88; fonttext2 = 28; colortext=[0.0 0.0 0.0];
    x0_colbar1 = 0.10; y0_colbar1 = 0.05; xl_colbar1 = 0.35; yl_colbar1 = 0.03;
    x0_colbar2 = 0.57;  y0_colbar2 = 0.05; xl_colbar2 = 0.35; yl_colbar2 = 0.03; 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %------------%%------------% Making figure %-----------%%------------%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %------------% Figure properties %------------%
    set(nfig,'position', [0, 100, 1000, 600]);
    
    p = panel(nfig);    
    p.pack({2/20 18/20}, 4);  
    p.fontsize=fontsize;
    
    
    %------------% Panel margins %------------%    

    p.de.margin=12;                                   %internal margins
    p(1,2).marginright = 20; p(2,2).marginright = 20; % Separate 2nd and 3rd columns
    p(1).marginbottom=8; p(2).margintop=0; 
    p.margin = [22 45 22 5];                        % margin[ = left bottom right top]

    
    %------------% Plot SL and SL rate 1st column %------------%
    p(1,1).select(); F01=gca;
    plot(F01,t,cos(par.omega*t),'--','linewidth',linew,'color',[0  0  0]); hold on; 
    plot(F01,t,sin(par.omega*t),':','linewidth',linew,'color',[0  0  0]); hold on; 
    plot(F01,t,(phip(:,end)-mean(phip(:,end)))./(max(phip(:,end))),'-','linewidth',linew,'color',[0 1 1]);
    set(F01,'Xlim',[trange(1) trange(end)],'Xtick',trange,'XTickLabel',[]);
    set(F01,'Ylim',[-1 1],'Ytick',[-1 0 1],'YTickLabel',{'$-1$' '0' '$1$'});
    set(F01,'Box','on');
    hold(F01,'on');
    grid(F01,'on');  

   
    %------------% Plot SL and SL rate 2nd column %------------%    
    p(1,2).select(); F01=gca;
    
    plot(F01,t,cos(par.omega*t),'--','linewidth',linew,'color',[0  0  0]); hold on; 
    plot(F01,t,sin(par.omega*t),':','linewidth',linew,'color',[0  0  0]); hold on; 
    plot(F01,t,(clp(:,end)-mean(clp(:,end)))./(max(clp(:,end))),'-','linewidth',linew,'color',[1 0 0]);
    set(F01,'Xlim',[trange(1) trange(end)],'Xtick',trange,'XTickLabel',[]);
    set(F01,'Ylim',[-1 1],'Ytick',[-1 0 1],'YTickLabel',[]);
    if (strcmp(par.RHSODE,'on')==1)
        set(F01,'YTickLabel',[]);
    else
        set(F01,'YTickLabel',[],'YAxisLocation','right');
    end
    set(F01,'Box','on');
    hold(F01,'on');
    grid(F01,'on');  

    %------------% Plot SL and SL rate 3rd column %------------%    
    p(1,3).select(); F01=gca;
    
    plot(F01,t,cos(par.omega*t),'--','linewidth',linew,'color',[0  0  0]); hold on; 
    plot(F01,t,sin(par.omega*t),':','linewidth',linew,'color',[0  0  0]); hold on; 
    plot(F01,t,(fp(:,end)-mean(fp(:,end)))./(max(fp(:,end))),'-','linewidth',linew,'color',[0 0 1]);
    set(F01,'Xlim',[trange(1) trange(end)],'Xtick',trange,'XTickLabel',[]);
    set(F01,'Ylim',[-1 1],'Ytick',[-1 0 1],'YTickLabel',[]);
    if (strcmp(par.RHSODE,'on')==1)
        set(F01,'YTickLabel',[]);
    else
        set(F01,'YTickLabel',[],'YAxisLocation','right');
    end
    set(F01,'Box','on');
    hold(F01,'on');
    grid(F01,'on');  

    
    %------------% Plot SL and SL rate 4th column %------------%    
    p(1,4).select(); F01=gca;
    
    plot(F01,t,cos(par.omega*t),'--','linewidth',linew,'color',[0  0  0]); hold on; 
    plot(F01,t,sin(par.omega*t),':','linewidth',linew,'color',[0  0  0]); hold on; 
    plot(F01,t,(fcp(:,end)-mean(fcp(:,end)))./(max(fcp(:,end))),'-','linewidth',linew,'color',[0.5 0 0.5]);
    set(F01,'Xlim',[trange(1) trange(end)],'Xtick',trange,'XTickLabel',[]);
    set(F01,'Ylim',[-1 1],'Ytick',[-1 0 1],'YTickLabel',{'$-1$' '0' '$1$'},'YAxisLocation','right');
    set(F01,'Box','on');
    hold(F01,'on');
    grid(F01,'on'); 

    
    %------------%%------------% Phip COLUMN %------------%%------------%
    phip0=phip;
    phip(abs(phip)<10^phiprange2(1))=10^phiprange2(1)*sign(phip(abs(phip)<10^phiprange2(1)));
    phip(abs(phip)>10^phiprange2(2))=10^phiprange2(2)*sign(phip(abs(phip)>10^phiprange2(2)));
    
    % Defining bounds
    logscale.min_sat  = phiprange2(1);
    logscale.max_sat  = phiprange2(2); 
    logscale.min      = min(min(log10(abs(phip0))));
    logscale.max      = max(max(log10(abs(phip0))));
    logscale.diff_sat = logscale.max_sat - logscale.min_sat + 1; %% +1 is to take into account lower bound for saturation
    logscale.end_sat  = round(logscale.min_sat - logscale.diff_sat);

    % Rescaling array
    phip_m = zeros(size(phip));
    phip_m(phip>0)  = log10(phip(phip>0));
    phip_m(phip<0)  = logscale.end_sat  - log10(-phip(phip<0)) + logscale.max_sat;
    phip_m(phip==0) = logscale.min_sat; 
    phip_m(phip_m<logscale.end_sat)=logscale.end_sat;

    ticks_rescaled =  logscale.end_sat:1:phiprange2(2);    

    str_colorbar_label = sprintf('$\\phi^{\\prime}$');      %% Label scaling 

    %------------% Plot Phip %------------%
    p(2,1).select(); F01=gca;   
    
    [FC,h]=contourf(T,Z,phip_m(1:end,:)',linspace(ticks_rescaled(1),ticks_rescaled(end),2*nc+1)); h.LineStyle='none'; hold on;
    caxis(F01,[min(ticks_rescaled) max(ticks_rescaled)]); 
    colormap(F01,cmap);
            
    set(F01,'Box','on');
    set(F01,'Xlim',[trange(1) trange(end)],'Xtick',trange,'XTickLabel',{'0' '$\frac{t_p}{2}$' '$t_p$' '$\frac{3 t_p}{2}$' '$2t_p$'});
    F01.XGrid = 'on'; F01.YGrid = 'on'; F01.Layer='top';
    text(-0.3,0.5,'$z$','Units','normalized','fontsize',24);
    text(xtext2,ytext2,str_colorbar_label,'Units','normalized','fontsize',fonttext2,'color',colortext); 
        
    % Managing colorbar
    colbar=colorbar(F01,'Location','southoutside','Position',[x0_colbar1 y0_colbar1 xl_colbar1 yl_colbar1],'Units','normalized');
    colbar.Limits=caxis;
    colbar.YAxisLocation = 'top';     
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
    colbar.Label.Interpreter='LaTex'; colbar.Label.String=str_colorbar_label; 
    colbar.Label.Units='normalized';  colbar.Label.Position=[-0.12 1.5]; % position relative to colorbar
    colbar.FontSize=15; colbar.Label.FontSize=18;
    
    
    %------------%%------------% clp COLUMN %------------%%------------%
    
    clp  = clp;
    clp0 = clp;
    clp(abs(clp)<10^clprange2(1))=10^clprange2(1)*sign(clp(abs(clp)<10^clprange2(1)));
    clp(abs(clp)>10^clprange2(2))=10^clprange2(2)*sign(clp(abs(clp)>10^clprange2(2)));
    
    % Defining bounds
    logscale.min_sat = clprange2(1);
    logscale.max_sat = clprange2(2); 
    logscale.min     = min(min(log10(abs(clp0))));
    logscale.max     = max(max(log10(abs(clp0))));
    logscale.diff_sat= logscale.max_sat - logscale.min_sat + 1; %% +1 is to take into account lower bound for saturation
    logscale.end_sat = round(logscale.min_sat - logscale.diff_sat);
    
    % Rescaling array
    clp_m = zeros(size(clp));
    clp_m(clp>0)  = log10(clp(clp>0));
    clp_m(clp<0)  = logscale.end_sat - log10(-clp(clp<0)) + logscale.max_sat ;
    clp_m(clp==0) = logscale.min_sat; 

    ticks_rescaled =  logscale.end_sat:1:clprange2(2);    

    str_colorbar_label = sprintf('$c_l^{\\prime}$'); %% Label scaling 
    
%     %----------%  Plot clp  %----------%
     p(2,2).select(); F01=gca;
    
    [FC,h]=contourf(T,Z,clp_m(1:end,:)',linspace(ticks_rescaled(1),ticks_rescaled(end),2*nc+1)); h.LineStyle='none'; hold on;
    caxis(F01,[min(ticks_rescaled) max(ticks_rescaled)]); 
    colormap(F01,cmap);
    
    set(F01,'Box','on');
    set(F01,'Xlim',[trange(1) trange(end)],'Xtick',trange,'XTickLabel',{'0' '$\frac{t_p}{2}$' '$t_p$' '$\frac{3 t_p}{2}$' '$2t_p$'});
    set(F01,'YTickLabel',[]);
    set(F01,'YAxisLocation','right');
    F01.XGrid = 'on'; F01.YGrid = 'on'; F01.Layer='top';
    text(xtext2,ytext2,str_colorbar_label,'Units','normalized','fontsize',fonttext2,'color',colortext); 

    
    % Managing colorbar
    colbar=colorbar(F01,'Location','southoutside','Position',[x0_colbar1 y0_colbar1 xl_colbar1 yl_colbar1],'Units','normalized');
    colbar.Limits=caxis;
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
%     ielim = [2,4]; %% Tick Labels to eliminate
%         for i = 1:length(ielim); colbar.TickLabels{ielim(i)}=' '; colbar.TickLabels{length(colbar.TickLabels)-ielim(i)+1}=' '; end;
    colbar.Label.Interpreter='LaTex'; colbar.Label.String=str_colorbar_label; 
    colbar.Label.Units='normalized';  colbar.Label.Position=[-0.12 -0.2]; % position relative to colorbar
    colbar.FontSize=15; colbar.Label.FontSize=18;

    
    %------------%%------------% Fp COLUMN %------------%%------------%
    fp0=fp;
    fp(abs(fp)<10^fprange2(1))=10^fprange2(1)*sign(fp(abs(fp)<10^fprange2(1)));
    fp(abs(fp)>10^fprange2(2))=10^fprange2(2)*sign(fp(abs(fp)>10^fprange2(2)));
    
    % Defining bounds
    logscale.min_sat  = fprange2(1);
    logscale.max_sat  = fprange2(2); 
    logscale.min      = min(min(log10(abs(fp0))));
    logscale.max      = max(max(log10(abs(fp0))));
    logscale.diff_sat = logscale.max_sat - logscale.min_sat + 1; %% +1 is to take into account lower bound for saturation
    logscale.end_sat  = round(logscale.min_sat - logscale.diff_sat);

    % Rescaling array
    fp_m = zeros(size(fp));
    fp_m(fp>0)  = log10(fp(fp>0));
    fp_m(fp<0)  = logscale.end_sat  - log10(-fp(fp<0)) + logscale.max_sat;
    fp_m(fp==0) = logscale.min_sat; 
    fp_m(fp_m<logscale.end_sat)=logscale.end_sat;

    ticks_rescaled =  logscale.end_sat:1:fprange2(2);    

    str_colorbar_label = sprintf('$f^{\\prime}$');      %% Label scaling 


     %----------%  Plot fp  %----------%
     p(2,3).select(); F01=gca;

    [FC,h]=contourf(T,Z,fp_m(1:end,:)',linspace(ticks_rescaled(1),ticks_rescaled(end),2*nc+1)); h.LineStyle='none'; hold on;
    caxis(F01,[min(ticks_rescaled) max(ticks_rescaled)]); 
    colormap(F01,cmap);

    set(F01,'Box','on');
    set(F01,'Xlim',[trange(1) trange(end)],'Xtick',trange,'XTickLabel',{'0' '$\frac{t_p}{2}$' '$t_p$' '$\frac{3 t_p}{2}$' '$2t_p$'});
    set(F01,'YTickLabel',[]);
    F01.XGrid = 'on'; F01.YGrid = 'on'; F01.Layer='top';
    text(xtext2,ytext2,str_colorbar_label,'Units','normalized','fontsize',fonttext2,'color',colortext); 

    % Managing colorbar
    colbar=colorbar(F01,'Location','southoutside','Position',[x0_colbar2 y0_colbar2 xl_colbar2 yl_colbar2],'Units','normalized');
    colbar.Limits=caxis;
    colbar.Ticks = ticks_rescaled;
    colbar.YAxisLocation = 'top';     
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
    colbar.Label.Interpreter='LaTex'; colbar.Label.String=str_colorbar_label; 
    colbar.Label.Units='normalized';  colbar.Label.Position=[-0.12 1.5]; % position relative to colorbar
    colbar.FontSize=15; colbar.Label.FontSize=18;
    %ielim


    %------------%%------------% Fcp COLUMN %------------%%------------%
    fcp0=fcp;
    fcp(abs(fcp)<10^fcprange2(1))=10^fcprange2(1)*sign(fcp(abs(fcp)<10^fcprange2(1)));
    fcp(abs(fcp)>10^fcprange2(2))=10^fcprange2(2)*sign(fcp(abs(fcp)>10^fcprange2(2)));

    % Defining bounds
    logscale.min_sat  = fcprange2(1);
    logscale.max_sat  = fcprange2(2); 
    logscale.min      = min(min(log10(abs(fcp0))));
    logscale.max      = max(max(log10(abs(fcp0))));
    logscale.diff_sat = logscale.max_sat - logscale.min_sat + 1; %% +1 is to take into account lower bound for saturation
    logscale.end_sat  = round(logscale.min_sat - logscale.diff_sat);

    % Rescaling array
    fcp_m = zeros(size(fcp));
    fcp_m(fcp>0)  = log10(fcp(fcp>0));
    fcp_m(fcp<0)  = logscale.end_sat  - log10(-fcp(fcp<0)) + logscale.max_sat;
    fcp_m(fcp==0) = logscale.min_sat; 
    fcp_m(fcp_m<logscale.end_sat)=logscale.end_sat;

    ticks_rescaled =  logscale.end_sat:1:fcprange2(2);    

    str_colorbar_label = sprintf('$f_c^{\\prime}$');      %% Label scaling


     %----------%  Plot fcp  %----------%
     p(2,4).select(); F01=gca;

    [FC,h]=contourf(T,Z,fcp_m(1:end,:)',linspace(ticks_rescaled(1),ticks_rescaled(end),2*nc+1)); h.LineStyle='none'; hold on;
    caxis(F01,[min(ticks_rescaled) max(ticks_rescaled)]); 
    colormap(F01,cmap);

    set(F01,'Box','on');
    text(1.2,0.5,'$z$','Units','normalized','fontsize',24);
    set(F01,'Xlim',[trange(1) trange(end)],'Xtick',trange,'XTickLabel',{'0' '$\frac{t_p}{2}$' '$t_p$' '$\frac{3 t_p}{2}$' '$2t_p$'});
    set(F01,'YAxisLocation','right');
    F01.XGrid = 'on'; F01.YGrid = 'on'; F01.Layer='top';
    text(xtext2,ytext2,str_colorbar_label,'Units','normalized','fontsize',fonttext2,'color',colortext); 


    % Managing colorbar
    colbar=colorbar(F01,'Location','southoutside','Position',[x0_colbar2 y0_colbar2 xl_colbar2 yl_colbar2],'Units','normalized');
    colbar.Limits=caxis;
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
    colbar.Label.Interpreter='LaTex'; colbar.Label.String=str_colorbar_label; 
    colbar.Label.Units='normalized';  colbar.Label.Position=[-0.12 -0.2]; % position relative to colorbar
    colbar.FontSize=15; colbar.Label.FontSize=18;
        
    
end