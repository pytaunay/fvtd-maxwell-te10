%%% Function flux_matrix
%%% Calculates the matrix from the split flux approach
%%% The matrix depends on the normal vector to a surface and the medium

function A = flux_matrix(n,mu,eps)

    nx = n(1,1);
    ny = n(1,1);
    nz = n(1,1);

    c = 1/sqrt(mu*eps); 

    A = zeros(6,6);

    % Upper left block (also lower right block)
    UL = [(ny^2+nz^2), -nx*ny, -nx*nz;
        -nx*ny, (nx^2+nz^2), -nz*ny;
        -nx*nz, -nz*ny, (nx^2+ny^2)];
    UL = UL*c;

    % Upper right block
    UR = [0,nz/eps,-ny/eps;
        -nz/eps,0,nx/eps;
        ny/eps,-nx/eps,0];

    % Lower left
    LL = [0,nz/mu,-ny/mu;
        -nz/mu,0,nx/mu; 
        ny/mu,-nx/mu,0];

    % Fill top left, bottom right
    A = kron(Ab,eye(2));
    
    % Fill top right and botom left
    A(1:3,4:6) = UR;
    A(4:6,1:3) = LL;
    

    %%% Result
    A = 1/2*A;

end
