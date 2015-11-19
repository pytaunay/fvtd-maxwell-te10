%%% Function eh_pec
%%% 
%%% Applies perfect electric boundary conditions to get the flux at a given
%%% finite volume surface

function ehs = eh_pec(Ei,Hi,Y,cpmat,idx)

    if idx == 1 || idx == 2
        dir = 1;
    elseif idx == 3 || idx == 4
        dir  = 2;
    elseif idx==5 || idx==6;
        dir =3;
    end

    Es = zeros(3,1);
    Hs = zeros(3,1);

    Es(dir,1) = Ei(dir,1);

    Htemp = Y*cpmat{dir}*cpmat{dir}*Ei+cpmat{dir}*Hi;

    switch dir
        case 1
            Hs(2,1) = Htemp(3,1);
            Hs(3,1) = -Htemp(2,1);
        case 2
            Hs(2,1) = -Htemp(3,1);
            Hs(3,1) = Htemp(2,1);
        case 3
            Hs(1,1) = -Htemp(3,1);
            Hs(3,1) = Htemp(1,1);
        case 4
            Hs(1,1) = Htemp(3,1);
            Hs(3,1) = -Htemp(1,1);
        case 5
            Hs(1,1) = Htemp(2,1);
            Hs(2,1) = -Htemp(1,1);
        case 6
            Hs(1,1) = -Htemp(2,1);
            Hs(2,1) = Htemp(1,1);
    end

    ehs = cat(1,Es,Hs);

end
    
