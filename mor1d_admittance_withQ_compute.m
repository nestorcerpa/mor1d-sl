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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----------%    Compute admittances    %----------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%----------% Parameters for computation of admittance %----------%

%%%         Wet          
Q_array  = [1e4     2e4     4e4     6e4     8e4     1e5     2e5     4e5     6e5     8e5     1e6]; 
nQs = length(Q_array);
tp_array = [1:1:50     52:2:100     105:5:150 ]; ntps = length(tp_array);

%----------% Loop over models %----------%
par.verb = "off";
for iQ = 1:nQs
    
    fprintf('\n ### Calculating admittance for model series : %2d ... \n', iQ);
    par.Q = Q_array(iQ); 

    fprintf(' --->  Q = %6.1e \n',par.Q)
    
    for itp = 1:ntps
        
        fprintf('%3d/%3d  ',itp,ntps);
        
        %----------% Update forcing period %----------%
        par.tp = tp_array(itp)*1e3; % [yr]
        par    = get_dimensionless_parameters(par);  % Updating dimensionless parameters (for omega);
        zarray = linspace(0,1,par.nz); 
         
        %----------%%----------%%----------%%----------%
        %----------% Calculating solution   %----------%
        %----------%%----------%%----------%%----------%

        %----------% Calculate steady state %----------%
        [MFields.cs,MFields.phi,~]         = mean_analytical(zarray,par); 
        %----------% Get other steady-state variables %----------%
        [MFields.W,MFields.q,MFields.qc] = get_other_mfields(MFields.phi,MFields.cs,par);
        
        %----------% Calculate fluctuations %----------%
        [FFields.csh,FFields.phih,~]       = fluctuations(zarray,par);
        %----------% Get other fluctuating variables %----------%
        [FFields.Wh,FFields.qh,FFields.qch] = get_other_ffields(MFields.phi,MFields.cs,FFields.phih,FFields.csh,par); 

        %----------%%----------%%----------%%----------%
        %----------%    Saving results      %----------%
        %----------%%----------%%----------%%----------%

        % Parameters array
        data.par_array(iQ,itp)       = par;
        % Fields at top of the column
        data.MFieldsTop.phi(iQ,itp)  = MFields.phi(end); 
        data.MFieldsTop.cs(iQ,itp)   = MFields.cs(end);
        data.MFieldsTop.q(iQ,itp)    = MFields.q(end);
        data.MFieldsTop.qc(iQ,itp)   = MFields.qc(end);
        data.FFieldsTop.phih(iQ,itp) = FFields.phih(end); 
        data.FFieldsTop.csh(iQ,itp)  = FFields.csh(end);
        data.FFieldsTop.qh(iQ,itp)   = FFields.qh(end);
        data.FFieldsTop.qch(iQ,itp)  = FFields.qch(end);      
        
    end
    
end
data.Q_array  = Q_array; 
data.tp_array = tp_array;
fprintf("\n\n ... DONE\n\n");

%----------% Saving data %----------%
save('mor1d_admittance_withQ.mat','data');
