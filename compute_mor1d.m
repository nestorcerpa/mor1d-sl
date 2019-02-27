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

%----------% Add paths to functions %----------%
restoredefaultpath;
addpath(genpath([pwd,'/src/']),genpath([pwd,'/external-functions/'])); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----------% Initialize parameters %----------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

par = input_parameters();

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
[~,~,MFields.cs,MFields.phi,~,~] = mean_analytical(zarray,par); %calcuate base state analytically

%----------% Get other mean variables %----------%
[MFields.W,MFields.f,MFields.fc] = get_other_mfields(MFields.phi,MFields.cs,par);
    
%----------% Plot mean state %----------%
nfig=nfig+1; figure(nfig); 
plot_meanfields(nfig,MFields,zarray,par,linew,fontsize);

%----------% Calculate fluctuations %----------%
[~,~,FFields.csh,FFields.phih] = fluctuations(zarray,par);
[FFields.Wh,FFields.fh,FFields.fch] = get_other_ffields(MFields.phi,MFields.cs,FFields.phih,FFields.csh,par); 

%----------% Plot perturbed state %----------%
nfig=nfig+1; figure(nfig); 
plot_fluctuationsh(nfig,MFields,FFields,zarray,par,linew,fontsize);

nfig=nfig+1; figure(nfig);
plot_fluctuations_ztspace(nfig,MFields,FFields,par,linew,fontsize,nperiods);

