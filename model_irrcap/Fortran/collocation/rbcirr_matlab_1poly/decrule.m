function [kdec cdec]=decrule(kgrid, theta, A_sol, order, kstart, kend)

global BETA ALPHA DELTAK GAMMAC PSI PHI iss


nk = length(kgrid);
ntheta = length(theta);

cdec=zeros(nk,ntheta);
kdec=zeros(nk,ntheta);

for this_k_pos =1:nk
    for this_theta_pos = 1:ntheta
        
        a_state = [A_sol(:,this_theta_pos)]';
        efunc = a_state*chebypol(kgrid(this_k_pos),order,kstart,kend);
        
        
        cdec(this_k_pos,this_theta_pos) = efunc^(-1/GAMMAC);
        kdec(this_k_pos,this_theta_pos) = invcfunc(theta(this_theta_pos),kgrid(this_k_pos),cdec(this_k_pos,this_theta_pos));
        
        % next verify guess for lambdak
        if kdec(this_k_pos,this_theta_pos)-(1-DELTAK)*kgrid(this_k_pos)<PHI*iss
            kdec(this_k_pos,this_theta_pos) = (1-DELTAK)*kgrid(this_k_pos)+PHI*iss;
            cdec(this_k_pos,this_theta_pos) = cfunc(theta(this_theta_pos),kgrid(this_k_pos), kdec(this_k_pos,this_theta_pos));
            
        end
    end
end

kdec = transpose(kdec);
cdec = transpose(cdec);