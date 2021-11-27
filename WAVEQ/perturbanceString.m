close all
clearvars
clc

% Perturbance over a string
L = 1; % meters
Nx = 1000;
c = 340;
dt = L/Nx/c;
C = 1;
T = 0.003;

x = 0:L/Nx:L;
sigma = 0.02;
mu = L/2;
I = exp(-(x-mu).^2/(2*sigma^2));

solverAndViz(I, 0, 0, c, L, dt, C, T, 10)