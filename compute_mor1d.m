close all; clear all;

%----------% Add paths to functions %----------%
restoredefaultpath;
addpath([pwd,'/src/'],genpath([pwd,'/external-functions/'])); 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----------% Initialize parameters %----------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

par = initialize_parameters();

%----------% Physical parameters %----------%
par.U0 = 5.0;            % Spreading-rate [cm/yr]
par.W0 = par.U0*par.cmyr_to_ms*2.0/pi;

%----------% Fluid dynamics parameters %----------%
par.n = 3;              % Permeability-porosity exponent
par.k = 1.e-6;          % Permeability [m^2]

%----------% Carbon parameters %----------%
par.cs0_vol  = 1e-4;    % Initial concentration of volatile element in solid 
par.D_vol    = 1e-4;    % Partition coefficient of volatile element in solid
par.M_vol    = 1e6;     % Solidus changes due to changes in composition (value Crowley et al. 2015)


%----------% Spatial and time parameters  %----------%
par.H = -(par.TsP0 - par.T_pot - par.M_vol*par.cs0_vol)/(par.nu^-1*par.rhom*par.grav)*1e-3; % Column height == depth at which T=Tsol for a given TsP0 
par.t0   = par.H*1e3/par.W0;  % [s] 

par.nz  = 2000;  zarray=linspace(0,1,par.nz);
par.tp  = 100e3;    % Forcing period [yr] 

%----------% Get dimensionless parameters %----------%
par = get_nondim_parameters(par); 

%----------% Plotting parameters %----------%
nfig=0;
linew=3;
fontsize=18;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----------%    Solving problem     %----------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%----------% Calculate mean state %----------%
[~,~,MFields.c,MFields.phi,~,~] = mean_analytical(zarray,par); %calcuate base state analytically

%----------% Plot mean state %----------%
nfig=nfig+1; figure(nfig); 
plot_meanfields(nfig,MFields,zarray,par,linew,fontsize)

%----------% Calculate fluctuations %----------%
[~,~,FFields.ch,FFields.phih] = fluctuations(zarray,par);

%----------% Plot perturbed state %----------%
nfig=nfig+1; figure(nfig) 
plot_fluctuationsh(nfig,MFields,FFields,zarray,par,linew,fontsize);

