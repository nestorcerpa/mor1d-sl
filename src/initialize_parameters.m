function params=initialize_parameters() 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %----------%  Initializing parameters for 1DColumn models %----------%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %----------% Generic parameters %----------%
    params.grav = 9.81;                 % Gravity [m/s^-2]
    params.cmyr_to_ms = 0.316881e-9;
    params.sec_per_year = 31557600;
    
    %----------% Physical parameters
    params.U0   = 5.0;          % Spreading-rate [cm/yr]
    params.rhom = 3300 ;        % Mantle density [kg/m3]
    params.Drho = 500 ;         % Mantle-Melt density contrast [kg/m3]
    params.rhow = 1000 ;        % Fresh water density [kg/m3]
    
    %----------% Fluid dynamics parameters
    params.n    = 2;          % Permeability-porosity exponent
    params.k    = 1e-6;       % Permeability C*d^2 [m^2]
    params.mu   = 1.0;        % Fluid viscosity [Pa s]
    
    %----------% Thermodynamic parameters %----------%
    params.T_pot = 1623.0;      % Mantle potential temperature [K]
    params.T_0   = 273;         % Ocean floor temperature      [K]
    params.TsP0  = 1503;        % solidus at P=0 [K]

    params.kappa = 1e-6;        % Thermal diffusivity
    params.nu  = 1./(60e-9);    % Clausius-Clapeyron slope of mantle (value Crowley et al. 2015)
    params.cp  = 1000;          % Specific heat capacity
    params.Lat = 6e5 ;          % Latent Heat

    %----------% Carbon parameters 
    params.cs0_vol  = 1e-4;     % Initial concentration of volatile element in solid 
    params.D_vol    = 1e-5;     % Partition coefficient of volatile element in solid
    params.M_vol    = 1.0 ;     % Solidus changes due to changes in composition (value Crowley et al. 2015)
    
    %----------% Spatial and temporal parameters               
    params.tp    = 40e3; % Period of SL changes [yr]
    params.ntime = 2000; % Time-steps
    
    params.nz      = 1000; % Number of points in z_array
    
    %----------% Calculating other reference parameters 
    params.H    = 100;
    params.W0   = params.U0*params.cmyr_to_ms*2/pi;
    params.t0   = params.H*1e3/params.W0;  % [s] 
    
    
    %----------% Miscellaneous options %----------%
    params.verb='on';               % Turning on/off displays
    params.ztype = 'linspace';      % Type of points distribution in z-array 
    params.RHSODE='on';             % Option to turn off RHS in ODEs solver

    %----------% ODE solver options %----------%
    params.ODEopts.RelTol = 1e-10;
    params.ODEopts.AbsTol = 1e-12;

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %----------%  Initializing various matlab options %----------%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
    %----------% Latex font %----------%
    if (strcmp(get(groot,'defaulttextinterpreter'),'latex')~=1); set(groot,'defaulttextinterpreter','latex'); end;
    if (strcmp(get(groot,'defaultAxesTickLabelInterpreter'),'latex')~=1); set(groot,'defaultAxesTickLabelInterpreter','latex'); end;
    if (strcmp(get(groot,'defaultLegendInterpreter'),'latex')~=1); set(groot,'defaultLegendInterpreter','latex'); end;
 

end
