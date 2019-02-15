function [xh,yh,ch,phih] = fluctuations(z_out,par)

%----------%%----------%%----------%%----------%%----------%
%PERTURBATION Calculates base state variables and derivative
%   Inputs
%       z_out is set of depths at which outputs are evaluated (can be a scalar or vector)
%       D is effective partition coefficient D=D/(1-D)
%       G=$\Gamma^*$
%       M=$\mathcal{M}$
%       Q=$\mathcal{Q}$
%       w=$\omega$ (non-dimensional frequency)
%   Outputs
%       xh=$\hat{x}$=M/c;
%       yh=$\hat{y}$
%       phih=$\hat{phi}$
%----------%%----------%%----------%%----------%%----------%

D = par.Deff; 
G = par.G;
M = par.M;
n = par.n;
Q = par.Q;
w = par.omega;

%----------% specifiy ode45 options (relatively tight tolerances)
options = odeset('AbsTol',par.ODEopts.AbsTol,'RelTol',par.ODEopts.RelTol,'Stats',par.verb);

%----------% calculate base state derivatives (needed for boundary conditions)
z=0;
[~,~,~,~,dx0,dy0] = mean_analytical(z,par);

%----------% calculate perturbed state
sol = ode45(@ode,[0 1],[-dx0;-dy0],options);

%----------% output
xh=deval(z_out,sol,1);
yh=deval(z_out,sol,2);
[~,~,~,phi,~,~] = mean_analytical(z_out,par);
dQdphi =  Q*(n*phi.^(n-1).*(1-phi).^2 - 2*(1-phi).*phi.^n); 
phih=yh./(1+dQdphi);
ch = xh/M; 

    function deriv = ode(z,y)
        xh = y(1);
        yh = y(2);
        [x,y,~,phi,dx,dy] = mean_analytical(z,par);
        dQdphi =  Q*(n*phi^(n-1)*(1-phi)^2 - 2*(1-phi)*phi^n);
        phih=yh/(1+dQdphi);
        dxh = (D+x+y)^(-1) *( 1i*w*G*x - xh*(1i*w*(D+phi+x)+dy) - yh*dx);
        dyh = 1i*w*(xh-G-phih) +dxh;
        deriv = [dxh;dyh];
    end

 
end

