close all
clearvars
clc

% Triangular wave (pulled string guitar)
L = 0.75; % meters
x0 = 0.2*L; 
a = 0.005;
freq = 440;
wavelength = 2*L;
c = freq*wavelength;
omega = 2*pi*freq;
num_periods = 2;
T = 2*pi/omega*num_periods;
Nx = 50;
% Choose dt the same as the stability limit for Nx=50
dt = L/Nx/c;
C = 1;

I = zeros(1,Nx+1);
x = 0:L/Nx:L;

for i = 1:Nx+1
    if (x(i) < x0)
        I(i) = a*x(i)/x0;
    else
        I(i) = a/(L-x0)*(L-x(i));
    end    
end

solverAndViz(I, 0, 0, c, L, dt, C, T, 1)