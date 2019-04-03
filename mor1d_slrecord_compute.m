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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----------% Reading input parameters %----------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

par = input_parameters();

%----------% Other parameters %----------%
z_out = 1; % z-coordinate at which we compute the fluctuations

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----------% Extract SL record within time-window (time) %----------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

time_window = [-1000. 0.];  % Time window for extracting data in SL-record
tp_min      = 1;            % Minimum sampling period in Fourier Transform

[SL] = ReadSLrecord('Siddall_2010.txt',time_window,tp_min);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----------% Computation of response to SL record %----------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%----------% Calculate steady state %----------%
[MFields.cs,MFields.phi,~] = mean_analytical(z_out,par); 
%----------% Get other steady-state variables %----------%
[MFields.W,MFields.q,MFields.qc] = get_other_mfields(MFields.phi,MFields.cs,par);

nperiod = length(SL.dfs.period); 

MODEL.phih_top  = zeros(nperiod,1); 
MODEL.clh_top   = zeros(nperiod,1); 
MODEL.qh_top    = zeros(nperiod,1); 
MODEL.qch_top   = zeros(nperiod,1); 

%----------% Calculate sampled response %----------%
fprintf('\n \t Starting calculation ... \n');   
nprog = 1;     par.verb  = 'off';
for iperiod = 1:nperiod

    xprog = iperiod/nperiod*100.;
    if (xprog >= nprog*5.0)
        fprintf('\t %4.1f..%',xprog); drawnow;
        nprog = nprog + 1;
    end

    %----------% Updating period %----------%
    par.tp    = SL.dfs.period(iperiod)*1e3; % [yr]
    par.omega = 2.*pi/par.tp*par.t0;

    %----------% Calculate fluctuating variables %----------%
    [FFields.csh,FFields.phih,~] = fluctuations(z_out,par);   % calculate fluctuating fields at z_out;
    %----------% Get other fluctuating variables %----------%
    [FFields.Wh,FFields.qh,FFields.qch] = get_other_ffields(MFields.phi,MFields.cs,FFields.phih,FFields.csh,par); 

    %----------% Saving fields at top of the column %----------%
    MODEL.phih_top(iperiod,1) = FFields.phih(end);
    MODEL.clh_top(iperiod,1)  = FFields.csh(end)/par.D_vol;
    MODEL.qh_top(iperiod,1)   = FFields.qh(end);
    MODEL.qch_top(iperiod,1)  = FFields.qch(end);

end
fprintf('\n \t ... End calculation \n');    

MODEL.SL        = SL;
MODEL.MFields   = MFields;
MODEL.par       = par;

datafile = 'mor1d_slrecord.mat';
save(datafile,'MODEL','-v7.3');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%     Other Functions     %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [SL] = ReadSLrecord(filename,timew,tp_min);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% READSLRECORD Reads sea-level record from .txt file 
%   Inputs
%      filename : name of the .txt file containing the sl record 
%      timew    : time window on which extract data from the sl record
%      tp_min   : minimum period in Fourier Transform
%   Outputs
%       SL      : structure array containing information on time-series sampling and 
%                 SL-Fourier Transform
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    fprintf('\n.... Reading Time Series of SL in %s...\n',filename);
    SLrecord_all = table2array(readtable(filename, 'HeaderLines', 1));

    % Extract Time-Series of SL between t1 and t2 (t1>0, t2>0 and t1<t2)
    [~,idx_start] = min(abs((SLrecord_all(:,1))-timew(1)));
    [~,idx_end] = min(abs((SLrecord_all(:,1))-timew(2)));
    
    SLrecord_window(:,1)  = SLrecord_all(idx_start:idx_end,1);
    SLrecord_window(:,2)  = SLrecord_all(idx_start:idx_end,2) - mean(SLrecord_all(idx_start:idx_end,2)); % SL relative to present-day SL

    % Extract arrays 
    time        = SLrecord_window(:,1);
    SLrecord    = SLrecord_window(:,2);

    % Mirror the time series
    SLrecord   = [flipud(SLrecord(2:end)); SLrecord];
    time = [time(1:end-1)+time(1); time];

    % Compute the fourier series and frequency vector
    SL.dfs = dft(SLrecord);
    k      = 1:length(SL.dfs.alpha);           % cycles
    SL.dfs.period = (time(end) - time(1))./k;  % 

    % Reconstruct SL discarding highest frequencies
    I            = find(SL.dfs.period >= tp_min);
    SL.dfs.alpha = SL.dfs.alpha(I);
    SL.dfs.beta  = SL.dfs.beta(I);
    SL.dfs.gamma = SL.dfs.gamma(I);
    SL.dfs.power = SL.dfs.power(I);
    SL.dfs.period= SL.dfs.period(I);
    SL.sl        = idft(SL.dfs);
    SL.time      = linspace(time(1),time(end),2*length(I)+1); 

end