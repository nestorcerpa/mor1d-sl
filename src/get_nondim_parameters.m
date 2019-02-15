function par=get_nondim_parameters(par)

    %----------%  Deff : effective partition coefficient D=D/(1-D)
    par.Deff    = par.D_vol/(1-par.D_vol);
    
    %----------%  G=$\Gamma^*$    
    par.G       = par.rhom*par.grav*par.H*1e3*par.cp/(par.nu*par.Lat);
    
    %----------%  M=$\mathcal{M}$
    par.M       = par.cp*abs(par.M_vol)*par.cs0_vol/par.Lat;
    
    %----------%  Q=$\mathcal{Q}$   
    par.Q       = par.Drho*par.grav*par.k/(par.mu*par.W0);

    if strcmp(par.verb,'on')==1
       fprintf('\n\n ###### Dimensionless parameters ##### \n')
       fprintf('  H = %4.1f, Deff = %6.4e, G = %6.2e, M = %6.2e, Q = %6.2e \n\n',par.H,par.Deff,par.G,par.M,par.Q); 
    end
  
    %----------%  \omega $   
    par.omega = 2*pi/(par.tp*par.sec_per_year)*par.t0;
    
end