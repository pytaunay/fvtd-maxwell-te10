%%% Function eh_pec
%%% 
%%% Applies perfect electric boundary conditions to get the flux at a given
%%% finite volume surface
%%%
%%% Flux along nk is 
%%% F(H).n = -n x H/eps0 = -1/eps0*(Y*n x n x E + n x H)
%%% F(E).n = n x E/mu = 0/mu for PEC


function ehs = eh_pec(Ei,Hi,Y,cpmat,idx)

    Es = zeros(3,1);

    Hs = Y*cpmat{idx}*cpmat{idx}*Ei+cpmat{idx}*Hi;


    ehs = cat(1,-Hs,Es);

end
    
