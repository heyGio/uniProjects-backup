function [y] = dattorro(x, Fs, depth, freq, blend, ff, fb, mod)
%DATTORRO  Industry standard effect structure
%
% J. Dattoro. Effect design, part 2: Delay-line modulation and chorus. 
% in J. Audio Eng. Soc., 45(10):764-788, October 1997

% length of the input signal
N = length(x);

% alloc the delay sequence
delay_seq = zeros(N, 1);

% compute the 'depth' parameter in (fractional) samples
depth_samples = depth * Fs;

% build the delay sequence
if(strcmp(mod,'sin'))
    
    % synthesize the modulating sinusoidal signal
    n = (1:N)/Fs;
	delay_seq = depth_samples + depth_samples * sin(2*pi*freq*n);
    
elseif(strcmp(mod,'noise'))
    
    % load the precomputed low-pass noise 
    ModLPN = load('MOD_lowpass_noise');
    field = strcat('lowpass_noise_', num2str(freq), 'Hz');
    lowpass_noise = ModLPN.(field);
    
    % synthesize the modulating noise signal
    delay_seq = depth_samples + depth_samples * lowpass_noise;
    
end

% alloc the output signal
y = zeros(N,1);
    
% delay line length
L = ceil(max(delay_seq))+3;
% alloc the delay line
delay_line = zeros(L,1);

for i=1:N
       
    % integer delay 
    d_int = floor(delay_seq(i));
    % fractionary delay
    d_frac = delay_seq(i) - d_int;
    % read a sample from the delay line ...
    %(using second-order FIR Lagrange interpolation filter with three coefficients)
    order = 2;
    h_int = zeros(order+1,1);
    for nn = 0:order
        k = (0:order)';
        k = k(k~=nn);
        h_int(nn+1) = prod((d_frac-k)./(nn-k));
    end

    x_d = h_int' * delay_line(d_int + 1: d_int + 3);

     
    
    % feedback the delayed input
    fb_comp = x_d * fb;
    
    % update the delay line
    delay_line(2:L) = delay_line(1:end-1); 
    delay_line(1) = x(i) - fb_comp;

    % feedforward component
    ff_comp = x_d * ff;
    
    % blend component
    blend_comp = delay_line(1) * blend;
    
    % generate the i-th sample of the output signal
    y(i) = ff_comp + blend_comp;
   
end

end