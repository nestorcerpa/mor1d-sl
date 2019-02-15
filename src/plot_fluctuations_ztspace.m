function plot_fluctuations_ztspace(nfig,MFields,FFields,zarray,par,linew,fontsize,nperiods);
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %------------%%------------% Defining arrays for z-t plots %------------%%------------%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    t = linspace(0,nperiods*par.tp*1e3*par.sec_per_year)./par.t0; % Non-dimensional time array
     
    phip  = real(FFields.phih.*exp(1i*par.omega*t));              % time-dependent fluctuating porosity
    clp   = real(FFields.ch./par.D_vol.*exp(1i*par.omega*t));     % time-dependent fluctuating liquid concentration
    
    slev   = real(exp(1i*par.omega*t));
    slrate = real(1i*par.omega*exp(1i*par.omega*t));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %------------%%------------% Making figure %-----------%%------------%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %------------% Initiale figure %------------%
    set(numfig,'position', [0, 100, 1000, 600]);
    
    p = panel(numfig);    
    p.pack({2/20 18/20}, 4);  
    p.fontsize=fontsize;
    
    
    %------------% Panel margins %------------%    

    p.de.margin=12;                                   %internal margins
    p(1,2).marginright = 20; p(2,2).marginright = 20; % Separate 2nd and 3rd columns
    p(1).marginbottom=8; p(2).margintop=0; 
    p.margin = [22 45 22 5];                        % margin[ = left bottom right top]

    
    %------------% Plot SL and SL rate 1st column %------------%
    p(1,1).select(); F01=gca;
    plot()
    
    

end