function plot_fluctuationsh(nfig,MFields,FFields,zarray,par,linew,fontsize)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plotting the amplitude and phase of \hat{c} and \hat{phi}
    %
    % *this function uses external-package panel-2.12
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
    plot(F01,log10(abs(FFields.ch)),zarray,'k','linewidth',linew);
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
    plot(F01,wrapTo2Pi(-angle(FFields.ch)),zarray,'k','linewidth',linew);
    xlabel(F01,'$\psi_{\hat{c}}$','Fontsize',fontsize);
    ylabel(F01,'$z$'); 
    set(F01,'xlim',[0 2*pi],'Xtick',[0 pi/2 pi 3*pi/2 2*pi],'Xticklabel',{'$0$','$\frac{\pi}{2}$','$\pi$','$\frac{3\pi}{2}$','$2\pi$'})
    set(F01,'ylim',[0 1],'Ytick',0:0.25:1);
    set(F01,'yaxislocation','right');
    set(F01,'Box','on');
    grid(F01,'on'); set(F01,'XMinorGrid','on','YMinorGrid','on');  
    
end