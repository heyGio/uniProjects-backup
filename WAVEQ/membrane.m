close all
clearvars
clc

Lx = 5;
Ly = 5;
dx = 0.1;
dy = dx;
nx = fix(Lx/dx);
ny = fix(Ly/dy);
x = linspace(0, Lx, nx);
y = linspace(0, Ly, ny);

T = 10;

wn = zeros(nx,ny);
wnm1 = wn;
wnp1 = wn;

C = 0.5;
c = 1;
dt = C*dx/c;

t = 0;

while (t < T) 
    wn(:,[1 end]) = 0;
    wn([1 end],:) = 0;
    
    t = t+dt;
    wnm1 = wn;
    wn = wnp1;
    
    %source
    %sourceX = 25; sourceY = 25;
    sourceX = 25; sourceY = 10;
    wn(sourceX,sourceY) = dt^2*20*sin(30*pi*t/20);
    
    for i=2:nx-1
        for j=2:ny-1
            wnp1(i,j) = 2*wn(i,j) - wnm1(i,j) ...
                        + C^2 *(wn(i+1,j) + wn(i,j+1) - 4*wn(i,j)... 
                        + wn(i-1,j) + wn(i,j-1));
        end
    end
    
    clf;
    subplot(2,1,1);
    imagesc(x,y,wn'); colorbar; caxis([-0.02 0.02]);
    
    
    subplot(2,1,2);
    mesh(x, y, wn'); colorbar; caxis([-0.02 0.02]);
    axis([0 Lx 0 Ly -0.05 0.05]);
    shg;
    pause(0.01);
    
end