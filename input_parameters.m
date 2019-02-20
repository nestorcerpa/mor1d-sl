function par=input_parameters() 

    %----------%%----------%%----------%%----------%
    %----------% ## INPUT PARAMETERS ## %----------%
    %----------%%----------%%----------%%----------%

    %----------% Physical parameters %----------%
    par.W0    = 2.0;    % Upwelling rate  [cm/yr]
    par.D_vol = 1e-4;   % Partition coefficient of carbon
    par.Fmax  = 0.2;    % Maximum degree of melting
    par.n     = 2;      % Exponent for permeability-porosity relationship
    par.RHSODE='on';             % Option to turn off RHS in ODEs solver

    %----------% Spatial and time parameters %----------%
    par.tp   = 100.e3;  % Period of forging             [yr]
    par.Hdry = 65.;     % Height of dry-melting column  [km]
    par.H    = 130.;    % Height of melting column      [km]

    %----------% Scaling parameters for fluctuations %----------%
    par.rhow = 1000;    % [kg/m3]
    par.rhom = 3300;    % [kg/m3]
    par.S0 = 50;        % [m]

    %----------% Input dimensionless parameters %----------%
    par.Q       = 1e5;

    %----------% ODE solver options %----------%
    par.ODEopts.RelTol = 1e-10;
    par.ODEopts.AbsTol = 1e-12;

    %----------% Miscellaneous options %----------%
    par.verb='on';               % Turning on/off displays
    
    %----------%%----------%%----------%%----------%
    %----------% ## OTHER PARAMETERS ## %----------%
    %----------%%----------%%----------%%----------%
    
    %----------% Derive other reference parameters %----------%
    par.t0   = par.H*1e3/(par.W0*1e-2); % time scale [yr]

    %----------% Define other dimensionless parameters %----------%
    par.Deff    = par.D_vol/(1-par.D_vol);
    par.G       = par.Fmax*(par.H/par.Hdry);         
    par.M       = par.Fmax*(par.H/par.Hdry - 1);
    par.omega   = 2.*pi/par.tp*par.t0;
    par.delta0  = par.S0/(par.H*1e3)*par.rhow/par.rhom;

    %----------% Spatial and time arrays %----------%
    par.nz = 2000;
    par.ntime = 1000;
 
    %----------% Other useful parameters %----------%
    par.cmyr_to_ms = 0.316881e-9;
    par.sec_per_year = 31557600;
        
    
    %----------%%----------%%----------%%----------%%----------%
    %----------%  ##  Modifying matlab options  ##  %----------%
    %----------%%----------%%----------%%----------%%----------%    
    
    %----------% Latex font %----------%
    if (strcmp(get(groot,'defaulttextinterpreter'),'latex')~=1); set(groot,'defaulttextinterpreter','latex'); end;
    if (strcmp(get(groot,'defaultAxesTickLabelInterpreter'),'latex')~=1); set(groot,'defaultAxesTickLabelInterpreter','latex'); end;
    if (strcmp(get(groot,'defaultLegendInterpreter'),'latex')~=1); set(groot,'defaultLegendInterpreter','latex'); end;

    
    %----------%%----------%%----------%%----------%%----------%
    %----------%  ##  Display input parameters  ##  %----------%
    %----------%%----------%%----------%%----------%%----------%
    if (strcmp(par.verb,'on')==1)
       fprintf('\n\n #-----# D = %4.1e ; Gamma = %4.2f ; M = %4.2f ; Q = %5.1e ; omega = %5.1f; #-----# \n\n',par.D_vol,par.G,par.M,par.Q,par.omega); 
    end
    
end


