function [x,y,cs,phi,dx,dy] = mean_analytical(z,par)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % MEAN_ANALYTICAL : Calculatesmean variables and derivatives
    %   Inputs 
    %       z is depth (can be a scalar or vector)
    %       par.Deff : effective partition coefficient D=D/(1-D)
    %       par.G=$\Gamma^*$
    %       par.M=$\mathcal{M}$
    %       par.Q=$\mathcal{Q}$
    %    Outputs 
    %       cs  : mean solid concentration
    %       phi : mean porosity
    %       x   : change of variable x = Mc
    %       y   : change of variable y = phi + q(phi)
    %       dx  : spatial-derivative of x
    %       dy  : spatial-derivative of y
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    D = par.Deff; 
    G = par.G;
    M = par.M;
    n = par.n; 
    Q = par.Q;  

    %----------% Calculate base state functions using quadratic formula 
    f=D+G*z-M;
    x=0.5*(-f+sqrt(f.^2+4*D*M));
    y=G*z-M+x;

    %----------% Determine phi (requires solution of an algebraic equation) 
    phi=zeros(size(y));
    for yt=y
        if n==2 %special case can be done with roots function
            pol_p = [Q, -2*Q, Q, 1,-yt];
            r_p=roots(pol_p);
            r_p = real(r_p(abs(imag(r_p))<1e-14));
            r_p=r_p(r_p<=1 & r_p>=0);
            if numel(r_p)~=1
                disp(r_p)
                error('incorrect number of roots found in the range (0,1)')
            end
            phi(yt==y)=r_p;
        elseif n==3 %special case can be done with roots function
            pol_p = [Q, -2*Q, Q, 0, 1,-yt];
            r_p=roots(pol_p);
            r_p = real(r_p(abs(imag(r_p))<1e-14));
            r_p=r_p(r_p<=1 & r_p>=0);
            if numel(r_p)~=1
                disp(r_p)
                error('incorrect number of roots found in the range (0,1)')
            end
            phi(yt==y)=r_p;
        else %fallback for n not integer
            phi(yt==y)=fzero(@(p) Q*p^n*(1-p)^2+p -yt,[0 1]);
        end
    end

    %----------% Calculate concentration (change variables)
    cs=x/M;

    %----------% Calculate derivatives
    dx=-x*G./(2*x+f); 
    dy=G+dx;

end
