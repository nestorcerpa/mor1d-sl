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

Q_array = [par.Q par.Q/4]; % Use two values of \mathcal{Q}
nQs = length(Q_array);

tp_array = [1:1:50     52:2:100     105:5:200]; % Define forcing periods
ntps = length(tp_array);

%----------% Loop over models %----------%
par.verb = "off";
for iQ = 1:nQs
    
    fprintf('\n ### Calculating admittance for Q=%4.1e ... \n', Q_array(iQ));
    par.Q = Q_array(iQ); 
    
    for itp = 1:ntps
        
        fprintf('%3d/%3d  ',itp,ntps);
        
        %----------% Update forcing period %----------%
        par.tp    = tp_array(itp)*1e3; % [yr]
        par.omega = 2.*pi/par.tp*par.t0;
         
        %----------%%----------%%----------%%----------%
        %----------% Calculating solution   %----------%
        %----------%%----------%%----------%%----------%

        %----------% Calculate mean variables %----------%
        [~,~,MFields.cs,MFields.phi,~,~]     = mean_analytical(zarray,par); 
        %----------% Get other mean variables %----------%
        [MFields.W,MFields.f,MFields.fc]    = get_other_mfields(MFields.phi,MFields.cs,par);
        %----------% Calculate fluctuating variables %----------%
        [~,~,FFields.csh,FFields.phih]       = fluctuations(zarray,par);
        %----------% Get other fluctuating variables %----------%
        [FFields.Wh,FFields.fh,FFields.fch] = get_other_ffields(MFields.phi,MFields.cs,FFields.phih,FFields.csh,par); 

        %----------%%----------%%----------%%----------%
        %----------%    Saving results      %----------%
        %----------%%----------%%----------%%----------%

        data.par_array(iQ,itp)       = par;
        data.MFieldsTop.phi(iQ,itp)  = MFields.phi(end); 
        data.MFieldsTop.cs(iQ,itp)   = MFields.cs(end);
        data.MFieldsTop.f(iQ,itp)    = MFields.f(end);
        data.MFieldsTop.fc(iQ,itp)   = MFields.fc(end);
        data.FFieldsTop.phih(iQ,itp) = FFields.phih(end); 
        data.FFieldsTop.csh(iQ,itp)  = FFields.csh(end);
        data.FFieldsTop.fh(iQ,itp)   = FFields.fh(end);
        data.FFieldsTop.fch(iQ,itp)  = FFields.fch(end);
        
    end
    
end
data.Q_array  = Q_array; 
data.tp_array = tp_array;
fprintf("\n\n ... DONE\n\n");

%----------% Saving data %----------%
save('column_admittance.mat','data');