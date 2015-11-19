%%% Author Pierre-Yves Taunay
%%% Nov. 2015
%%%
%%% Goal
%%% Simulate a straight section of WR112 waveguide
%%% using the FVTD method for Maxwell's equations in vac.
%%%
%%% Details
%%% The Maxwells equations are put in a conservative form, and the usual
%%% finite volume treatment is performed. 
%%% The numerical flux at the cell boundaries is calculated using a Steger
%%% and Warming flux splitting method. Time stepping is done with RK4 
%%% and explicit Euler

clear all;

%%% Define the geometry , initial power, constants
a = 28.4988e-3;
b = 12.6238e-3;
L = 5e-2;
f = 8e9; % Hz
P10 = 1; % Watts 
tmax = 1/f; % N time steps 

mu = 4*pi*1e-7;
eps = 8.851e-12;
c = 1/sqrt(mu*eps);

% Number of discretization points
Nx = 15;
Ny = 5;
Nz = 30;
Nt = 100; % Number of 
Ns = 6; % number of sides

dx = a/Nx;
dy = b/Ny;
dz = L/Nz;
dt = 1/(f*Nt); 

dV = dx*dy*dz;

% Derived values
om = 2*pi*f; 
beta = beta_te10(om,mu,eps,a);
Y = sqrt(eps/mu);

% Get coefficient for E and H
A10 = sqrt(4*pi^2*P10/(om*mu*a^3*b)*1/real(beta));

%%% Construct initial values for E and H
Ue = cell(Nx,Ny,Nz);
Uh = cell(Nx,Ny,Nz);
[Ue{:}] = deal(zeros(3,1));
[Uh{:}] = deal(zeros(3,1));

Uen = cell(Nx,Ny,Nz);
Uhn = cell(Nx,Ny,Nz);
[Uen{:}] = deal(zeros(3,1));
[Uhn{:}] = deal(zeros(3,1));

Ueall = cell(Nx,Ny,Nz,ceil(tmax/dt));
[Uen{:}] = deal(zeros(3,1));


%%% Create some temporary matrix holders

% Matrix of surface normals
Nmat = [-1,1,0,0,0,0;
        0,0,-1,1,0,0;
        0,0,0,0,-1,1];

% Matrix to extract perpendicular components
xpmat = zeros(6,6);
xpmat(2,2) = 1;
xpmat(3,3) = 1;
xpmat(5,5) = 1;
xpmat(6,6) = 1;

ypmat = zeros(6,6);
ypmat(1,1) = 1;
ypmat(3,3) = 1;
ypmat(4,4) = 1;
ypmat(6,6) = 1;

zpmat = zeros(6,6);
zpmat(1,1) = 1;
zpmat(2,2) = 1;
zpmat(4,4) = 1;
zpmat(5,5) = 1;

Pmat = cell(6,1);
[Pmat{1:2}] = deal(xpmat);
[Pmat{3:4}] = deal(ypmat);
[Pmat{5:6}] = deal(zpmat);

% Matrix to perform cross product with a unit vector
cpmat = cell(6,1);
xpmat = zeros(3,3);
xpmat(2,3) = -1;
xpmat(3,2) = 1;

ypmat = zeros(3,3);
ypmat(1,3) = 1;
ypmat(3,1) = -1;

zpmat = zeros(3,3);
zpmat(1,2) = -1;
zpmat(2,1) = 1;

cpmat{1} = xpmat;
cpmat{2} = -xpmat;
cpmat{3} = ypmat;
cpmat{4} = -ypmat;
cpmat{5} = zpmat;
cpmat{6} = -zpmat;


alpha_mat = eye(6);
alpha_mat(1,1) = 1/eps;
alpha_mat(2,2) = 1/eps;
alpha_mat(3,3) = 1/eps;

alpha_mat(4,4) = 1/mu;
alpha_mat(5,5) = 1/mu;
alpha_mat(6,6) = 1/mu;

%%% Time stepping
t = 0;
nt = 1;
tic()

%%% Original values
for i=1:Nx
    for j=1:Ny
        for k=1:Nz
            Ueall{i,j,k,nt} = Ue{i,j,k};
        end
    end
end

%while t < tmax
for m=1:5
	t = t+dt;
	%%% Calculate the flux
	for i=1:Nx
		for j=1:Ny
            for k=1:Nz
                Un = zeros(6,1);
    
                % Get my neighbors and cell surface
                NN_idx = cell(Ns,1);
                NN_idx{1} = [i-1,j,k,dy*dz]'; 
                NN_idx{2} = [i+1,j,k,dy*dz]'; 
                NN_idx{3} = [i,j-1,k,dx*dz]'; 
                NN_idx{4} = [i,j+1,k,dx*dz]'; 
                NN_idx{5} = [i,j,k-1,dx*dy]'; 
                NN_idx{6} = [i,j,k+1,dx*dy]'; 

                % Get my E-H field
                Ei = Ue{i,j,k}; 
                Hi = Uh{i,j,k}; 

                %%% Get the flux for each neighbor
                % Initialize flux vector
                F = zeros(6,1);
                for l=1:Ns
                    idx = NN_idx{l};
                    coord = idx(1:3);
                    S = idx(4);
                    
                    if coord(1) <= 0 || coord(1) >= Nx || coord(2) <= 0 || coord(2) >= Ny || coord(3) >= Nz
                        %%% Boundary condition: perfect electric conductor
                        EHs = eh_pec(Ei,Hi,Y,cpmat,l); 
                        %F = F+alpha_mat*eh_flux(EHs(1:3,1),EHs(4:6,1))*Nmat(:,l)*S;
                        F = F+alpha_mat*eh_flux(zeros(3,1),zeros(3,1))*Nmat(:,l)*S;
                    elseif coord(3) <= 0 %|| coord(3) >= Nz
                        %%% Boundary condition: TE10 mode excitation
                        % Center of surface coordinate
                        xs = dx*(i-1) + dx/2; 
                        if(coord(3) <= 0)
                            zs = 0;
                        else
                            zs = L;
                        end

                        EHs = eh_te10(xs,zs,t,a,A10,om,mu,eps,beta);

                        F = F+alpha_mat*eh_flux(EHs(1:3,1),EHs(4:6,1))*Nmat(:,l)*S;
                    else
                        %%% General case
                        % Get the neighbor's fields
                        El = Ue{coord(1),coord(2),coord(3)};
                        Hl = Uh{coord(1),coord(2),coord(3)};
                        
                        % From the definition of the flux function 
                        Fmaxw = 1/2*alpha_mat*eh_flux(El,Hl)*Nmat(:,l);
                        % Cross product term
                        Fcp = alpha_mat*Pmat{l}*[Ei-El;Hi-Hl];
                    
                        F = F+(Fmaxw+Fcp)*S;
                    end
                end

                %%% Euler to update the fields
                Un = cat(1,Ei,Hi) - 1/dV*dt*F; 

                Uen{i,j,k} = Un(1:3,1); 
                Uhn{i,j,k} = Un(4:6,1); 

            end
        end
    end
    Ue = Uen;
    Uh = Uhn;

    nt = nt+1;
    for i=1:Nx
        for j=1:Ny
            for k=1:Nz
                Ueall{i,j,k,nt} = Ue{i,j,k};
            end
        end
    end
end
    
toc()


