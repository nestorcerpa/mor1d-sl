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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----------% Reading input parameters %----------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

par = input_parameters();

zarray = linspace(0,1,par.nz);

%%%%---- Parameters for 2-d model approximation 

drho = 500;           % solid-liquid density contrast 
rhol = par.rhom-drho; % liquid density
rhoc = 2900.;         % crust density
par.alpha = 30.;      % angle of decompaction channel
par.xf    = 33.2;     % focusing distance xf
hf    = par.xf*tan(par.alpha*pi/180.); % depth decompating channel at x=xf
nz    = 30;
zf_array  = linspace(par.H-hf,par.H,nz)/par.H;
tan_alpha = tan(par.alpha*pi/180.);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----------%    Compute admittances    %----------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%----------% Parameters for computation of admittance %----------%

%%%         Wet         Wet         Wet         Dry         Dry Basal
Q_array  = [par.Q       4*par.Q     par.Q/4     par.Q       par.Q]; 
Harray   = [par.H       par.H       par.H       par.Hdry    par.Hdry];
Gammap   = {'on'        'on'        'on'        'on'        'off'}; 
nQs = length(Q_array);

tp_array = [1:1:50     52:2:100     105:5:200   300:50:500]; % Define forcing periods
ntps = length(tp_array);

%----------% Loop over models %----------%
par.verb = "off";
for iQ = 1:nQs
    
    fprintf('\n ### Calculating admittance for model series : %2d ... \n', iQ);
    par.Q      = Q_array(iQ); 
    par.Gammap = Gammap{iQ};
    par.H      = Harray(iQ);
    par=get_dimensionless_parameters(par);  % Updating dimensionless parameters 

    fprintf(' --->  H = %4.1f km; Q = %6.1e ; Gammap : %s ; t0 = %6.1e \n ',par.H,par.Q,par.Gammap,par.t0);
    
    %----------% Calculate mean variables %----------%
    [MFields.cs,MFields.phi,~]         = mean_analytical(zarray,par); %calculate base state analytically
    %----------% Get other mean variables %----------%
    [MFields.W,MFields.q,MFields.qc] = get_other_mfields(MFields.phi,MFields.cs,par);

    %----------% Calculate 2-d approx %----------%
    MFields.Rmor = 0.0; MFields.Rcmor = 0.0; MFields.Hcru = 0.0; 
    zf_array = linspace(par.H-hf,par.H,nz)/par.H;
    for i=2:length(zf_array)
        dz  = zf_array(i)-zf_array(i-1);
        [cs_i,phi_i,~] = mean_analytical(zf_array(i),par); %calculate base state analytically 
        [W_i,q_i,qc_i] = get_other_mfields(phi_i,cs_i,par);  % Get other variables    
        MFields.Rmor  = MFields.Rmor  + 2.0*q_i*dz/tan_alpha;
        MFields.Rcmor = MFields.Rcmor + 2.0*qc_i*dz/tan_alpha;
        MFields.Hcru  = MFields.Hcru  + 2.0*q_i/(pi/2)*rhol/rhoc*dz/tan_alpha;
        tau_array(i)  = 0.0; % tranport-time for melt originated at z = zf_array(i)
        for j = i:length(zf_array)
            dzj  = zf_array(j)-zf_array(j-1);
            [cs_j,phi_j,~] = mean_analytical(zf_array(j),par);  % Calculate base state analytically 
            w_j = par.Q*phi_j^(par.n-1)*(1-phi_j)^2;            % steady-state melt velocity
            tau_array(i) = tau_array(i) + 1./w_j*dzj;
        end
    end
    
    
    for itp = 1:ntps
        
        fprintf(' %3d/%3d  ',itp,ntps);
        
        %----------% Update forcing period %----------%
        par.tp    = tp_array(itp)*1e3; % [yr]
        par=get_dimensionless_parameters(par);  % Updating dimensionless parameters (for omega)
         
        %----------%%----------%%----------%%----------%
        %----------% Calculating solution   %----------%
        %----------%%----------%%----------%%----------%

        %----------% Calculate fluctuating variables %----------%
        [FFields.csh,FFields.phih,~]       = fluctuations(zarray,par);
        %----------% Get other fluctuating variables %----------%
        [FFields.Wh,FFields.qh,FFields.qch] = get_other_ffields(MFields.phi,MFields.cs,FFields.phih,FFields.csh,par); 

        %----------% 2-d focusing model %----------%
        fprintf('.');
        FFields.Rmorh0= 0.0; FFields.Rcmorh0= 0.0;
        FFields.Rmorh = 0.0; FFields.Rcmorh = 0.0; 
        for i=2:length(zf_array)
            dz = zf_array(i)-zf_array(i-1); 
            [cs_i,phi_i,~]    = mean_analytical(zf_array(i),par); % calculate steady state analytically 
            [ch_i,phih_i]     = fluctuations(zf_array(i),par);
            [~,qh_i,qch_i]    = get_other_ffields(phi_i,cs_i,phih_i,ch_i,par); 
            tau = tau_array(i); 
            FFields.Rmorh0  = FFields.Rmorh0  + 2.0*qh_i*dz/tan_alpha;
            FFields.Rcmorh0 = FFields.Rcmorh0 + 2.0*qch_i*dz/tan_alpha;
            FFields.Rmorh   = FFields.Rmorh   + 2.0*exp(-1i*par.omega*tau)*qh_i*dz/tan_alpha;
            FFields.Rcmorh  = FFields.Rcmorh  + 2.0*exp(-1i*par.omega*tau)*qch_i*dz/tan_alpha;
        end
        fprintf('.');
        
        %----------%%----------%%----------%%----------%
        %----------%    Saving results      %----------%
        %----------%%----------%%----------%%----------%

        %----------% Parameters array
        data.par_array(iQ,itp)       = par;
        %----------% Fields at top of the column
        data.MFieldsTop.phi(iQ,itp)  = MFields.phi(end); 
        data.MFieldsTop.cs(iQ,itp)   = MFields.cs(end);
        data.MFieldsTop.q(iQ,itp)    = MFields.q(end);
        data.MFieldsTop.qc(iQ,itp)   = MFields.qc(end);
        data.FFieldsTop.phih(iQ,itp) = FFields.phih(end); 
        data.FFieldsTop.csh(iQ,itp)  = FFields.csh(end);
        data.FFieldsTop.qh(iQ,itp)   = FFields.qh(end);
        data.FFieldsTop.qch(iQ,itp)  = FFields.qch(end);
        %----------% 2-d approximation 
        data.MFieldsMOR.Hcru(iQ,itp)    = MFields.Hcru;
        data.MFieldsMOR.Rmor(iQ,itp)    = MFields.Rmor;
        data.MFieldsMOR.Rcmor(iQ,itp)   = MFields.Rcmor;
        data.FFieldsMOR.Rmorh0(iQ,itp)  = FFields.Rmorh0;
        data.FFieldsMOR.Rcmorh0(iQ,itp) = FFields.Rcmorh0;
        data.FFieldsMOR.Rmorh(iQ,itp)   = FFields.Rmorh;
        data.FFieldsMOR.Rcmorh(iQ,itp)  = FFields.Rcmorh;
        %----------% Save Bottom-to-surface Lag 
        phase_phih  =  unwrap(angle(FFields.phih)); 
        phase_csh   =  unwrap(angle(FFields.csh) ); 
        phase_qh    =  unwrap(angle(FFields.qh)  ); 
        phase_qch   =  unwrap(angle(FFields.qch) );
        % 1-d model
        data.lagBCtoSurf.phi(iQ,itp) = -(phase_phih(end)-phase_phih(1))*tp_array(itp)/(2*pi);
        data.lagBCtoSurf.cs(iQ,itp)  = -(phase_csh(end)-phase_csh(1))*tp_array(itp)/(2*pi);
        data.lagBCtoSurf.q(iQ,itp)   = -(phase_qh(end)-phase_qh(1))*tp_array(itp)/(2*pi);
        data.lagBCtoSurf.qc(iQ,itp)  = -(phase_qch(end)-phase_qch(1))*tp_array(itp)/(2*pi);
        % 2-d model
        data.lagBCtoSurfMOR.Rmor0(iQ,itp)   = data.lagBCtoSurf.q(iQ,itp)  - min(  abs((angle(FFields.qh(end)) -angle(FFields.Rmorh0))*par.tp*1e-3/(2*pi)), abs((angle(FFields.qh(end))  - (angle(FFields.Rmorh0)+2*pi))*par.tp*1e-3/(2*pi))   );
        data.lagBCtoSurfMOR.Rcmor0(iQ,itp)  = data.lagBCtoSurf.qc(iQ,itp) - min(  abs((angle(FFields.qch(end))-angle(FFields.Rcmorh0))*par.tp*1e-3/(2*pi)),abs((angle(FFields.qch(end)) - (angle(FFields.Rcmorh0)+2*pi))*par.tp*1e-3/(2*pi))  );
        data.lagBCtoSurfMOR.Rmor(iQ,itp)    = data.lagBCtoSurf.q(iQ,itp)  - min(  abs((angle(FFields.qh(end)) -angle(FFields.Rmorh))*par.tp*1e-3/(2*pi)),  abs((angle(FFields.qh(end))  - (angle(FFields.Rmorh)+2*pi))*par.tp*1e-3/(2*pi))    );
        data.lagBCtoSurfMOR.Rcmor(iQ,itp)   = data.lagBCtoSurf.qc(iQ,itp) - min(  abs((angle(FFields.qch(end))-angle(FFields.Rcmorh))*par.tp*1e-3/(2*pi)), abs((angle(FFields.qch(end)) - (angle(FFields.Rcmorh)+2*pi))*par.tp*1e-3/(2*pi))   );
        
    end
    
end
data.Q_array  = Q_array; 
data.tp_array = tp_array;
fprintf("\n\n ... DONE\n\n");

%----------% Computation of admittance with parameters of Burley and Katz 2015 %----------%

fprintf('\n ### Calculating admittance for models with parameters of Burley and Katz 2015 \n');
par = input_parameters(0);

%%% Burley and Katz 2015 parameters 
cmyr_to_ms = 0.316881e-9;
U0_BK15 = 3.0; k0_BK15 = 1e-12; drho_BK15 = 500.; n_BK15=3; phi0_BK15 = 0.01;
W0_BK15 = U0_BK15*2.0/pi; % in m/s
w0_BK15 = drho_BK15*10.0*k0_BK15/phi0_BK15/cmyr_to_ms;
Q_BK15 = (w0_BK15/W0_BK15)^n_BK15 * par.Fmax^(1-n_BK15);


par.W0 = W0_BK15;
par.n  = n_BK15; 
par.Q  = Q_BK15;
par.H  = par.Hdry;
Gammap   = {'on' 'off'}; 
nQs = length(Gammap);
par=get_dimensionless_parameters(par);  % Updating dimensionless parameters 
    

par.verb = 'off';
for iQ = 1:nQs
    
    fprintf('\n ### Calculating Fields for BK15-Model series %2d ... \n',iQ);
    par.Gammap = Gammap{iQ};
    
    par=get_dimensionless_parameters(par);  % Updating dimensionless parameters 

    fprintf(' --->  Gammap : %s \n ',par.Gammap);
    
    %----------% Calculate steady state %----------%
    [MFields.cs,MFields.phi,~]         = mean_analytical(zarray,par); 

    %----------% Get other steady-state variables %----------%
    [MFields.W,MFields.q,MFields.qc] = get_other_mfields(MFields.phi,MFields.cs,par);

    %----------% Calculate 2-d approx %----------%
    MFields.Rmor = 0.0; MFields.Rcmor = 0.0; MFields.Hcru = 0.0; 
    zf_array = linspace(par.H-hf,par.H,nz)/par.H;
    for i=2:length(zf_array)
        dz  = zf_array(i)-zf_array(i-1);
        [cs_i,phi_i,~] = mean_analytical(zf_array(i),par); %calculate base state analytically 
        [W_i,q_i,qc_i] = get_other_mfields(phi_i,cs_i,par);  % Get other variables    
        MFields.Rmor  = MFields.Rmor  + 2.0*q_i*dz/tan_alpha;
        MFields.Rcmor = MFields.Rcmor + 2.0*qc_i*dz/tan_alpha;
        MFields.Hcru  = MFields.Hcru  + 2.0*q_i/(pi/2)*rhol/rhoc*dz/tan_alpha;
        tau_array(i)  = 0.0; % tranport-time for melt originated at z = zf_array(i)
        for j = i:length(zf_array)
            dzj  = zf_array(j)-zf_array(j-1);
            [cs_j,phi_j,~] = mean_analytical(zf_array(j),par);  % Calculate base state analytically 
            w_j = par.Q*phi_j^(par.n-1)*(1-phi_j)^2;            % steady-state melt velocity
            tau_array(i) = tau_array(i) + 1./w_j*dzj;
        end
    end
    
    for itp = 1:ntps

        fprintf(' %3d/%3d  ',itp,ntps);

        %----------% Update forcing period %----------%
        par.tp  = tp_array(itp)*1e3; % [yr]
        par     = get_dimensionless_parameters(par);  % Updating dimensionless parameters (for omega)

        %----------%%----------%%----------%%----------%
        %----------% Calculating solution   %----------%
        %----------%%----------%%----------%%----------%

        %----------% Calculate fluctuations %----------%
        [FFields.csh,FFields.phih,~]       = fluctuations(zarray,par);
        %----------% Get other fluctuating variables %----------%
        [FFields.Wh,FFields.qh,FFields.qch] = get_other_ffields(MFields.phi,MFields.cs,FFields.phih,FFields.csh,par); 

        %----------% 2-d focusing model %----------%
        fprintf('.');
        FFields.Rmorh0= 0.0; FFields.Rcmorh0= 0.0;
        FFields.Rmorh = 0.0; FFields.Rcmorh = 0.0; 
        for i=2:length(zf_array)
            dz = zf_array(i)-zf_array(i-1); 
            [cs_i,phi_i,~]    = mean_analytical(zf_array(i),par); % calculate steady state analytically 
            [ch_i,phih_i]     = fluctuations(zf_array(i),par);
            [~,qh_i,qch_i]    = get_other_ffields(phi_i,cs_i,phih_i,ch_i,par); 
            tau = tau_array(i); 
            FFields.Rmorh0  = FFields.Rmorh0  + 2.0*qh_i*dz/tan_alpha;
            FFields.Rcmorh0 = FFields.Rcmorh0 + 2.0*qch_i*dz/tan_alpha;
            FFields.Rmorh   = FFields.Rmorh   + 2.0*exp(-1i*par.omega*tau)*qh_i*dz/tan_alpha;
            FFields.Rcmorh  = FFields.Rcmorh  + 2.0*exp(-1i*par.omega*tau)*qch_i*dz/tan_alpha;
        end
        fprintf('.');

        %----------%%----------%%----------%%----------%
        %----------%    Saving results      %----------%
        %----------%%----------%%----------%%----------%

        % Parameters array
        data.par_array_BK15(iQ,itp)       = par;
        % Fields at top of the column
        data.MFieldsTop_BK15.phi(iQ,itp)  = MFields.phi(end); 
        data.MFieldsTop_BK15.cs(iQ,itp)   = MFields.cs(end);
        data.MFieldsTop_BK15.q(iQ,itp)    = MFields.q(end);
        data.MFieldsTop_BK15.qc(iQ,itp)   = MFields.qc(end);
        data.FFieldsTop_BK15.phih(iQ,itp) = FFields.phih(end); 
        data.FFieldsTop_BK15.csh(iQ,itp)  = FFields.csh(end);
        data.FFieldsTop_BK15.qh(iQ,itp)   = FFields.qh(end);
        data.FFieldsTop_BK15.qch(iQ,itp)  = FFields.qch(end);
        % 2-d approximation 
        data.MFieldsMOR_BK15.Hcru(iQ,itp)    = MFields.Hcru;
        data.MFieldsMOR_BK15.Rmor(iQ,itp)    = MFields.Rmor;
        data.MFieldsMOR_BK15.Rcmor(iQ,itp)   = MFields.Rcmor;
        data.FFieldsMOR_BK15.Rmorh0(iQ,itp)  = FFields.Rmorh0;
        data.FFieldsMOR_BK15.Rcmorh0(iQ,itp) = FFields.Rcmorh0;
        data.FFieldsMOR_BK15.Rmorh(iQ,itp)   = FFields.Rmorh;
        data.FFieldsMOR_BK15.Rcmorh(iQ,itp)  = FFields.Rcmorh;
        % Bottom to Surface lag
        phase_phih  =  unwrap(angle(FFields.phih)); 
        phase_csh   =  unwrap(angle(FFields.csh) ); 
        phase_qh    =  unwrap(angle(FFields.qh)  ); 
        phase_qch   =  unwrap(angle(FFields.qch) ); 
        data.lagBCtoSurf_BK15.phi(iQ,itp) = -(phase_phih(end)-phase_phih(1))*tp_array(itp)/(2*pi);
        data.lagBCtoSurf_BK15.cs(iQ,itp)  = -(phase_csh(end)-phase_csh(1))*tp_array(itp)/(2*pi);
        data.lagBCtoSurf_BK15.q(iQ,itp)   = -(phase_qh(end)-phase_qh(1))*tp_array(itp)/(2*pi);
        data.lagBCtoSurf_BK15.qc(iQ,itp)  = -(phase_qch(end)-phase_qch(1))*tp_array(itp)/(2*pi);
        data.lagBCtoSurfMOR_BK15.q(iQ,itp)     = -( (phase_qh(end)+angle(FFields.Rmorh0)) - phase_qh(1) )*tp_array(itp)/(2*pi);
        data.lagBCtoSurfMOR_BK15.qc(iQ,itp)    = -( (phase_qch(end)+angle(FFields.Rcmorh0)) - phase_qch(1) )*tp_array(itp)/(2*pi);
        
    end
end
fprintf("\n\n ... DONE\n\n");


%----------% Saving data %----------%
save('mor1d_admittance.mat','data');
