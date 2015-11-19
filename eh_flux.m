%%% Function
%%% eh_flux
%%% returns the flux from the maxwell's equation in a conservative form

function eh_f = eh_flux(E,H)

    F = zeros(6,3);

    F(1,1) = 0;
    F(2,1) = H(3,1);
    F(3,1) = -H(2,1);
    F(4,1) = 0;
    F(5,1) = -E(3,1);
    F(6,1) = E(2,1);

    
    F(1,2) = -H(3,1);
    F(2,2) = 0;
    F(3,2) = H(1,1);
    F(4,2) = E(3,1);
    F(5,2) = 0;
    F(6,2) = -E(1,1);

    F(1,3) = H(2,1);
    F(2,3) = -H(1,1);
    F(3,3) = 0;
    F(4,3) = -E(2,1);
    F(5,3) = E(1,1);
    F(6,3) = 0;

    eh_f = F;
end
