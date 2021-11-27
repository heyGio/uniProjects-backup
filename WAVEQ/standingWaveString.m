close all
clearvars
clc

% Perturbance over a string
L = 1; % meters
Nx = 1000;
c = 340;
dt = L/Nx/c;
C = 1;
T = 0.05;

x = 0:L/Nx:L;
t = 0:dt:T;
Nt = length(t);


%simulation
%mode / frameSkip
%mode = 1; frameSkip = 150;
mode = 3; frameSkip = 50;

f = zeros(Nx, Nt);
f(round(Nx/2),:)=  cos(mode*pi/(L)*c*t);


solverAndViz(0, 0, f, c, L, dt, C, T, frameSkip)