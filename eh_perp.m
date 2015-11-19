%%% wFunction eh_perp
%%% returns [E_perp, H_perp] respective to a direction
%%% i.e. projection in a plane

function eh_p  = eh_perp(E,H,dir,mu,eps)

    Ep = E;
    Hp = H;

    Ep(dir,1) = 0;
    Hp(dir,1) = 0;

    alpha_mat = eye(6);
    alpha_mat(1,1) = 1/eps;
    alpha_mat(2,2) = 1/eps;
    alpha_mat(3,3) = 1/eps;

    alpha_mat(4,4) = 1/mu;
    alpha_mat(5,5) = 1/mu;
    alpha_mat(6,6) = 1/mu;

    eh_p = alpha_mat * cat(1,Ep,Hp);

end
