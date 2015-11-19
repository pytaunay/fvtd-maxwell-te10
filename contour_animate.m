[XG,ZG] = meshgrid(dx/2:dx:a,dz/2:dz:L);

figure();


Econtour = zeros(size(XG,1),size(XG,2));
for k=1:Nz
    for i=1:Nx
        Econtour(k,i) = sum(abs(Ueall{i,3,k,1}));
    end
end


hCont=surf(XG,ZG,Econtour);
axis([0,a,0,L]);
pause;


for t=2:tmax/dt
    Econtour = zeros(size(XG,1),size(XG,2));
    for k=1:Nz
        for i=1:Nx
            Econtour(k,i) = sum(abs(Ueall{i,3,k,t}));
        end 
    end

    set(hCont,'ZData',Econtour);    
    
    pause;
end
    