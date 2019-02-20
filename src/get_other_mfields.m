function [W,f,fc] = get_other_mfields(phi,c,par)

    Q = par.Q; n = par.n;
    
    %----------% Mean solid velocity %----------% 
    W  = 1 - Q .* (phi.^n.*(1-phi)); % \bar{W}
    
    %----------% Mean melt flux %----------% 
    f  = 1 - (1-phi).*W;
    
    %----------% Mean chemical flux %----------% 
    fc = f.*c/par.D_vol;
    
end