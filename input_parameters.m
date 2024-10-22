% % % MIT Licence
% % % 
% % % Copyright (c) 2019
% % %     Nestor G. Cerpa       (University of Montpellier) [nestor.cerpa@gm.univ-montp2.fr]
% % %     David W. Rees Jones   (University of Oxford)      [david.reesjones@earth.ox.ac.uk]
% % %     Richard F. Katz       (University of Oxford)      [richard.katz@earth.ox.ac.uk] 
% % % 
% % % Permission is hereby granted, free of charge, to any person obtaining a copy
% % % of this software and associated documentation files (the "Software"), to deal
% % % in the Software without restriction, including without limitation the rights
% % % to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% % % copies of the Software, and to permit persons to whom the Software is
% % % furnished to do so, subject to the following conditions:
% % % 
% % % The above copyright notice and this permission notice shall be included in all
% % % copies or substantial portions of the Software.
% % % 
% % % THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% % % IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% % % FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% % % AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% % % LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% % % OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% % % SOFTWARE.

function par=input_parameters(varargin);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT_PARAMETERS 
%   Input 
%       varargin{1} : = 0 to not restart the default path 
%   Output
%       par : array with model parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    fprintf('\n READING PARAMETERS IN INPUT_PARAMETERS() ...\n');

    %----------%%----------%%----------%%----------%
    %----------% MODEL INPUT PARAMETERS %----------%
    %----------%%----------%%----------%%----------%

    %----------% Physical parameters %----------%
    par.W0     = 2.0;    % Upwelling rate  [cm/yr]
    par.D_vol  = 1e-4;   % Partition coefficient of carbon
    par.Fmax   = 0.2;    % Maximum degree of melting
    par.n      = 2;      % Exponent for permeability-porosity relationship
    par.Gammap ='on';   % Option to run a basal-flux model like Burley and Katz (2015)

    %----------% Spatial and time parameters %----------%
    par.tp   = 100.e3;  % Period of forging             [yr]
    par.Hdry = 65.;     % Height of dry-melting column  [km]
    par.H    = 130.;    % Height of melting column      [km]

    %----------% Scaling parameters for fluctuations %----------%
    par.rhow = 1000;    % Sea-water density                 [kg/m3]
    par.rhom = 3300;    % Mantle density                    [kg/m3]
    par.S0 = 50;        % Amplitude of sea-level vairiation [m]

    %----------% Parameters for pseudo-2-d models %----------%
    par.alpha = 30.;                % angle of decompaction channel
    par.xf    = 33.5;               % focusing distance xf
    par.hf    = par.xf*tan(par.alpha*pi/180.); % depth decompating channel at x=xf
    
    %----------% Parameters for crust production %----------%
    par.drho  = 500;                % solid-melt density contrast
    par.rhol  = par.rhom-par.drho;  % melt density
    par.rhoc  = 2900;               % oceanic crust density    
    
    
    %----------% Input dimensionless parameters %----------%
    par.Q       = 1e5;  % $\mathcal{Q}$
    
    
    %----------%%----------%%----------%
    %----------%   OPTIONS  %----------%
    %----------%%----------%%----------%
    
    %----------% ODE solver options %----------%
    par.ODEopts.RelTol = 1e-10;
    par.ODEopts.AbsTol = 1e-12;

    %----------% Miscellaneous options %----------%
    par.verb='on';               % Verbosity
   

    %----------% Spatial and time arrays %----------%
    par.dz = 0.25;    % cell size in km
    par.ntime = 100; % number of time steps

        
    %----------% Add paths to functions %----------%
    if (nargin==0 || varargin{1}~=false);
        fprintf('...restoring default matlab path...')
        restoredefaultpath; addpath(genpath([pwd,'/src/']),genpath([pwd,'/external-functions/'])); 
         %----------% Derive dimensionless parameters %----------%
        par=get_dimensionless_parameters(par);
    end;
    
    %----------%%----------%%----------%%----------%%----------%
    %----------%  ##  Modifying matlab options  ##  %----------%
    %----------%%----------%%----------%%----------%%----------%    
    
    %----------% Latex font %----------%
    if (strcmp(get(groot,'defaulttextinterpreter'),'latex')~=1); set(groot,'defaulttextinterpreter','latex'); end;
    if (strcmp(get(groot,'defaultAxesTickLabelInterpreter'),'latex')~=1); set(groot,'defaultAxesTickLabelInterpreter','latex'); end;
    if (strcmp(get(groot,'defaultLegendInterpreter'),'latex')~=1); set(groot,'defaultLegendInterpreter','latex'); end;
        
end


