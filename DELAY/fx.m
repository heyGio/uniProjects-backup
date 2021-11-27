%------------------------------------------%
%        *** SSSP - HOMEWORK #4 ***        %
%------------------------------------------%
%   Audio effect based on a time-varying   %
%          fractional delay line           %
%------------------------------------------%
% Giovanni Affatato, Roberto Alessandri    %
%------------------------------------------%

clear; close all; clc;

%% Fx type
% - 'vibrato'
% - 'flanger'
% - 'chorus'
% - 'white_chorus'
% - 'doubling'
% - 'echo'
fx_type = 'echo' %#ok<NOPTS>

%% Read the input file
[x, Fs] = audioread('NylonGt1-C4-mono.wav');

%% FX parameters 
switch lower(fx_type)
    
    case {'vibrato'}
        depth=0.002; freq=2; blend=0; ff=1; fb=0; mod='sin';
            
    case {'flanger'}
        depth=0.005; freq=0.5; blend=0.707; ff=0.707; fb=-0.707; mod='sin';
        
    case {'chorus'}
        depth=0.015; freq=10; blend=1.; ff=0.707; fb=0; mod='noise';
               
    case {'white_chorus'}
        depth=0.015; freq=10; blend=0.707; ff=1.; fb=0.707; mod='noise';
        
    case {'doubling'}
        depth=0.07; freq=3; blend=0.707; ff=0.707; fb=0.; mod='noise';
        
    case {'echo'}
        depth=0.4; freq=0; blend=1.; ff=0.8; fb=0.4; mod='sin';
    
    otherwise
        error('fx_type \"%s\" not found.', fx_type)
                
end

%% Apply FX
y = dattorro(x, Fs, depth, freq, blend, ff, fb, mod);

%% Avoid any (possible) clipping
y = rescale(y,-1.,1.);

%% Playback
%audiowrite([fx_type,'.wav'], y, Fs);
soundsc(y, Fs)

%% Read the reference audio file
dir_name = 'FX_ref';
addpath(dir_name);
[y_ref, ~] = audioread(fullfile(dir_name, strcat(fx_type,'.wav')));

%% Display the MSE
MSE = mean(abs(y-y_ref).^2);
MSE_str = sprintf('MSE: %g', MSE);
disp(MSE_str)
