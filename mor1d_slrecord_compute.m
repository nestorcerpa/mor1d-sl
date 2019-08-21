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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----------% Extract SL record within time-window (time) %----------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

time_window = [-1000. 0.];  % Time window for extracting data in SL-record
tp_min      = 1;            % Minimum sampling period in Fourier Transform

[SL] = ReadSLrecord('Siddall_2010.txt',time_window,tp_min);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----------% Computation of response to SL record %----------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

w0_ref    = par.Fmax^((par.n-1)/par.n) * par.Q^(1./par.n)*par.W0; % reference maximum melt velocity
w0_array  = [w0_ref 100. 200. 400. 600. 800. 1000.];              % [in cm/yr]

fprintf('\n\n'); par.verb = 'off';
for iw = 1:length(w0_array)

    fprintf('\n\n \t #------##------# CALCULATING SL-RECORD MODEL : %3d/%3d ...', iw,length(w0_array));
    fprintf(' ... with w0 = %5.1f cm/y',w0_array(iw))
    par.Q  = (w0_array(iw)/par.W0)^(par.n) * par.Fmax^(1-par.n);   % Updating \mathcal{Q} with new w_0
    par    = get_dimensionless_parameters(par);                    % Updating dimensionless parameters    
    zarray = linspace(0,1,par.nz);           
    tan_alpha = tan(par.alpha*pi/180.);
    
    %----------% Calculate steady state %----------%
    [MFields.cs,MFields.phi,~] = mean_analytical(zarray,par); 
    %----------% Get other steady-state variables %----------%
    [MFields.W,MFields.q,MFields.qc] = get_other_mfields(MFields.phi,MFields.cs,par);
    %----------% Calculate steady-state pseudo-2-d model %----------%
    [~,idx_zf] = min(abs(zarray-(1-par.hf/par.H)));
    MFields.Rmor = 0.0; MFields.Rcmor = 0.0; MFields.Hcru = 0.0; MFields.tau_array(1:idx_zf) = nan;
    for i=idx_zf+1:length(zarray)
        dz  = zarray(i)-zarray(i-1);
        MFields.Rmor  = MFields.Rmor  + 2*MFields.q(i)*dz/tan_alpha; % Factor 2 is for taking into account both sides of ridge's axis
        MFields.Rcmor = MFields.Rcmor + 2*MFields.qc(i)*dz/tan_alpha;
        MFields.Hcru  = MFields.Hcru  + 2*MFields.q(i)/(pi/2)*par.rhol/par.rhoc*dz/tan_alpha;
        MFields.tau_array(i)  = 0.0; % tranport-time for melt originated at z = zf_array(i)
        for j = i:length(zarray)
            dzj  = zarray(j)-zarray(j-1);
            [cs_j,phi_j,~] = mean_analytical(zarray(j),par);  % Calculate base state analytically 
            w_j = par.Q*phi_j^(par.n-1)*(1-phi_j)^2;          % steady-state melt velocity
            MFields.tau_array(i) = MFields.tau_array(i) + 1./w_j*dzj;
        end
    end
    fprintf("\n\t\t ... Calculated mean crustal thickness = %5.2f km \n",MFields.Hcru*par.H);
        
    nperiod = length(SL.dfs.period); 
    DATA.model{iw}.phih_top  = zeros(nperiod,1); 
    DATA.model{iw}.clh_top   = zeros(nperiod,1); 
    DATA.model{iw}.qh_top    = zeros(nperiod,1); 
    DATA.model{iw}.qch_top   = zeros(nperiod,1); 
    DATA.model{iw}.Rmorh0    = zeros(nperiod,1); DATA.model{iw}.Rcmorh0    = zeros(nperiod,1);
    DATA.model{iw}.Rmorh1    = zeros(nperiod,1); DATA.model{iw}.Rcmorh1    = zeros(nperiod,1);
    
    %----------% Calculate sampled response %----------%
    fprintf('\n \t Starting loop over periods ... \n');   
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
        [FFields.csh,FFields.phih,~] = fluctuations(zarray,par);   % calculate fluctuating fields at z_out;
        %----------% Get other fluctuating variables %----------%
        [FFields.Wh,FFields.qh,FFields.qch] = get_other_ffields(MFields.phi,MFields.cs,FFields.phih,FFields.csh,par); 
        %----------% Calculate fluctuating pseudo-2-d model %----------%
        Rmorh0 = 0.0; Rmorh1 = 0.0; Rcmorh0 = 0.0; Rcmorh1 = 0.0;
        for i=idx_zf+1:length(zarray)
            dz  = zarray(i)-zarray(i-1);
            qh_i = FFields.qh(i);
            qch_i = FFields.qch(i);
            % Instantaneous focusing     
            Rmorh0  = Rmorh0  + 2.0*qh_i*dz/tan_alpha;
            Rcmorh0 = Rcmorh0 + 2.0*qch_i*dz/tan_alpha;
            % Focusing time equal to steady-state melt transport time
            Rmorh1   = Rmorh1  + 2.0*exp(-1i*par.omega*MFields.tau_array(i))*qh_i*dz/tan_alpha;
            Rcmorh1  = Rcmorh1 + 2.0*exp(-1i*par.omega*MFields.tau_array(i))*qch_i*dz/tan_alpha;
        end
        
        %----------% Saving fields at top of the column %----------%
        DATA.model{iw}.MFields             = MFields;
        DATA.model{iw}.phih_top(iperiod,1) = FFields.phih(end);
        DATA.model{iw}.clh_top(iperiod,1)  = FFields.csh(end)/par.D_vol;
        DATA.model{iw}.qh_top(iperiod,1)   = FFields.qh(end);
        DATA.model{iw}.qch_top(iperiod,1)  = FFields.qch(end);
        DATA.model{iw}.Rmorh0(iperiod,1)   = Rmorh0;
        DATA.model{iw}.Rmorh1(iperiod,1)   = Rmorh1;
        DATA.model{iw}.Rcmorh0(iperiod,1)  = Rcmorh0;
        DATA.model{iw}.Rcmorh1(iperiod,1)  = Rcmorh1;
    end
    DATA.model{iw}.par = par;  
end
DATA.SL      = SL;
DATA.tp_min  = tp_min;
fprintf('\n \t ... DONE \n');  
fprintf('\n \t ... Saving results in mor1d_slrecord.mat \n'); 

datafile = 'mor1d_slrecord.mat';
save(datafile,'DATA','-v7.3');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%     Other Functions     %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function SL = ReadSLrecord(filename,timew,tp_min);

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