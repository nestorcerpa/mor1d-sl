function plot_meanfields(nfig,MFields,zarray,par,linew,fontsize)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plotting mean variables \bar{c},\bar{phi},\bar{f} and \bar{f_c}
    %
    % *this function uses external-package panel-2.12
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
    semilogx(F01,MFields.c,zarray,'k','linewidth',linew); hold(F01,'on');
    xlabel(F01,'$\overline{c}$','Fontsize',fontsize);
    set(F01,'ylim',[0 1],'Ytick',0:0.25:1)
    set(F01,'YtickLabel',[]);       
    set(F01,'Box','on');
    grid(F01,'on'); set(F01,'XMinorGrid','on','YMinorGrid','on'); 
    
    
    %----------% Melt flux %----------%
    f = 1. + (1.0 - MFields.phi).*(par.Q .* (MFields.phi.^par.n) .* (1.-MFields.phi) - 1.0);
    p(1,3).select(); F01=gca; 
    semilogx(F01,f,zarray,'k','linewidth',linew); hold(F01,'on');
    xlabel(F01,'$\overline{f}$','Fontsize',fontsize);
    set(F01,'ylim',[0 1],'Ytick',0:0.25:1)
    set(F01,'YtickLabel',[]);       
    set(F01,'Box','on');
    grid(F01,'on'); set(F01,'XMinorGrid','on','YMinorGrid','on');  
    
    %----------% Chemical flux %----------%
    fc = MFields.c/par.D_vol .* f;
    p(1,4).select(); F01=gca; 
    semilogx(F01,fc,zarray,'k','linewidth',linew); hold(F01,'on');
    xlabel(F01,'$\overline{f_c}$','Fontsize',fontsize);
    set(F01,'ylim',[0 1],'Ytick',0:0.25:1)
    set(F01,'yaxislocation','right');       
    set(F01,'Box','on');
    grid(F01,'on'); set(F01,'XMinorGrid','on','YMinorGrid','on');    
 
end