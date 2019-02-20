close all; clear all;

%----------% Add paths to functions %----------%
restoredefaultpath;
addpath([pwd,'/src/'],genpath([pwd,'/external-functions/'])); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----------% Initialize parameters %----------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

par = input_parameters();

%----------% Modify some parameters %----------%
par.tp = 100.e3; 

%----------% Spatial and time arrays %----------%
zarray = linspace(0,1,par.nz);

nperiods = 2; tf = nperiods*par.tp/par.t0;
time   = linspace(0,2*tf,par.ntime);

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
[FFields.Wh,FFields.fh,FFields.fch] = get_other_ffields(MFields.phi,MFields.c,FFields.phih,FFields.ch,par); 

%----------% Plot perturbed state %----------%
nfig=nfig+1; figure(nfig); 
plot_fluctuationsh(nfig,MFields,FFields,zarray,par,linew,fontsize);

nfig=nfig+1; figure(nfig);
plot_fluctuations_ztspace(nfig,MFields,FFields,par,linew,fontsize,nperiods);

