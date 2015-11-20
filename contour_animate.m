[XG,ZG] = meshgrid(dx/2:dx:a,dz/2:dz:L);

figure();


Econtour = zeros(size(XG,1),size(XG,2));
Hcontour = zeros(size(XG,1),size(XG,2));
for k=1:Nz
    for i=1:Nx
        Econtour(k,i) = sum(abs(Ueall{i,3,k,1}));
        Hcontour(k,i) = sum(abs(Uhall{i,3,k,1}));
    end
end

subplot(2,1,1);
eCont=surf(XG,ZG,Econtour);
axis([0,a,0,L]);
subplot(2,1,2);
hCont=surf(XG,ZG,Hcontour);
axis([0,a,0,L]);
pause;


for t=2:tmax/dt
    Econtour = zeros(size(XG,1),size(XG,2));
    Hcontour = zeros(size(XG,1),size(XG,2));
    for k=1:Nz
        for i=1:Nx
            Econtour(k,i) = sum(abs(Ueall{i,3,k,t}));
            Hcontour(k,i) = sum(abs(Uhall{i,3,k,t}));
        end 
    end

    set(eCont,'ZData',Econtour);    
    set(hCont,'ZData',Hcontour);   
    pause;
end
    