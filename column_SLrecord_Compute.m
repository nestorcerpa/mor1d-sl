close all; clear all;

%----------% Add paths to functions %----------%
restoredefaultpath;
addpath([pwd,'/src/'],genpath([pwd,'/external-functions/'])); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----------% Reading input parameters %----------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

par = input_parameters();

%----------% Other parameters %----------%
z_out = 1; % z-coordinate for computation of the results

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----------% Extract SL record within time-window (time) %----------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

time_window = [-1000. 0.];  % Time window for extracting data in SL-record
tp_min      = 1;            % Minimum sampling period in Fourier Transform

[SL,SLrecord_window] = Func_ReadSLrecord_cutF('Siddall_2010.txt',time_window,tp_min);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----------% Computation of response to SL record %----------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%----------% Calculate mean state %----------%
[~,~,MFields.cs,MFields.phi,~,~] = mean_analytical(z_out,par); % calculate mean state analytically
%----------% Get other mean variables %----------%
[MFields.W,MFields.f,MFields.fc] = get_other_mfields(MFields.phi,MFields.cs,par);

nperiod = length(SL.dfs.period); 

MODEL.phih_top  = zeros(nperiod,1); 
MODEL.clh_top   = zeros(nperiod,1); 
MODEL.fh_top    = zeros(nperiod,1); 
MODEL.fch_top   = zeros(nperiod,1); 

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
    [~,~,FFields.csh,FFields.phih] = fluctuations(z_out,par); % calculate fluctuating fields at z_out;
    %----------% Get other fluctuating variables %----------%
    [FFields.Wh,FFields.fh,FFields.fch] = get_other_ffields(MFields.phi,MFields.cs,FFields.phih,FFields.csh,par); 

    %----------% Saving fields at top of the column %----------%
    MODEL.phih_top(iperiod,1) = FFields.phih(end);
    MODEL.clh_top(iperiod,1)  = FFields.csh(end)/par.D_vol;
    MODEL.fh_top(iperiod,1)   = FFields.fh(end);
    MODEL.fch_top(iperiod,1)  = FFields.fch(end);

end
fprintf('\n \t ... End calculation \n');    

MODEL.SL        = SL;
MODEL.MFields   = MFields;
MODEL.par       = par;

datafile = 'column_SLRecord.mat';
save(datafile,'MODEL','-v7.3');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%     Other Functions     %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [SL,SLrecord_window] = Func_ReadSLrecord_cutF(filename,time_window,tp_min);

    fprintf('\n.... Reading Time Series of SL in %s...\n',filename);
    SLrecord_all = load(filename);

    % Extract Time-Series of SL between t1 and t2 (t1>0, t2>0 and t1<t2)
    [~,idx_start] = min(abs((SLrecord_all(:,1))-time_window(1)));
    [~,idx_end] = min(abs((SLrecord_all(:,1))-time_window(2)));
    
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
    SL.dfs.period = (time(end) - time(1))./k;  % kiloyears

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