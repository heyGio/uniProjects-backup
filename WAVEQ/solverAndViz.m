function [] = solverAndViz(I_u,V_u,f_u,c,L,dt,C,T, frame_skip)

Nt = round(T/dt);
t = 0:dt:T;
dx = dt*c/C;
Nx = round(L/dx);
x = 0:dx:L;
C2 = C^2; 

if f_u == 0
    f = zeros(Nx+1, Nt+1);
else
    f = f_u;
end



if I_u == 0
    I = zeros(1, Nx+1);
else
    I = I_u;
end

if V_u == 0
    V = zeros(1, Nx+1);
else
    V = V_u;
end

ris_mesh = zeros(Nt+2, Nx+1);

u_n = I;
ris_mesh(1, :) = u_n;


n = 1;
u = zeros(1,Nx+1);
for i = 2: Nx
    u(i) = u_n(i) - dt*V(i)  + 0.5*C2*(u_n(i-1) - 2*u_n(i) + u_n(i+1)) + 0.5*dt^2*f(i, n);
end
u(1) = 0; 
u(Nx+1) = 0;

u_nm1 = u_n;
u_n = u;
ris_mesh(2, :) = u_n;

for n = 2:Nt
    % Update all inner points at time t[n+1]
    for i = 2: Nx
        u(i) = - u_nm1(i) + 2*u_n(i) + C2*(u_n(i-1) - 2*u_n(i) + u_n(i+1)) + dt^2*f(i, n);
    end
    % Insert boundary conditions
    u(1) = 0; 
    u(Nx+1) = 0;
    
    
    u_nm1 = u_n;
    u_n = u;
    ris_mesh(n+1, :) = u_n;
end


figure
for i=1:frame_skip:Nt
    clf
    hold on
    
    
    plot(x, ris_mesh(i,:), 'LineWidth', 2)
    grid on

    axis([0 max(x) -max(max(ris_mesh)) max(max(ris_mesh))])
    drawnow
end

% figure;
% imagesc(ris_mesh);

end

