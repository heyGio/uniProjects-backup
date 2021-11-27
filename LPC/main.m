%----------------------------------------------%
%          *** DAAP - HOMEWORK #1 ***          %
%----------------------------------------------%
%           Linear Predictive Coding           %
%----------------------------------------------%
% Giovanni Affatato, Roberto Alessandri        %
%----------------------------------------------%
clear; close all; clc;

% order of the LPC filters
p = 64;

% control string for which of the three algorithms to apply; the variable 
% takes values in {"levinson", "steepest", "lms"}
algorithm = "levinson";

% read "speech.wav" 
[speech, Fs] = audioread("speech.wav");

% audio playback
input_player = audioplayer(speech, Fs);
% play(input_player); 

% define the analysis window length
duration = 40*10^-3;
dur_samples = duration * Fs;
n = ceil(log2(dur_samples));
M = 2^n;

% define the analysis window
window = hann(M);

% define the hop-size
Hop = M/2;

% define the minimum allowed Kronecker comb's period
Tkc_min = 1/200;

% initialize output signal to a zero vector
out = zeros(length(speech), 1);

% define total number of analysis windows
n_windows = floor((length(speech)-M)/Hop) + 1;


for i = 0:n_windows-1
    % window the speech signal
    frame = window .* speech(i*Hop + 1:i*Hop + M);
   
        
    % subtract the sample mean from the short-time speech segment
    sample_mean = mean(frame);
    frame = frame - sample_mean;  
        
    
    % compute the zero-crossing rate of the speech segment    
    zcr = (1/(M-1))*sum(abs(diff(sign(frame))));
    
    
    if algorithm == "levinson"
        [a, e] = levinson_durbin(frame, p);      
    elseif algorithm == "steepest"  
        [a, e] = steepest_descent(frame, p, 2000);
    elseif algorithm == "lms"
        [a, e] = lms_algorithm(frame, p, 0.03);
    end
    
    % compute the power of the prediction error 
    sigma_e2 = (1/M)*sum(e.^2);    
    
    if zcr > 0.2
       % compute the excitation signal
       white_noise = randn(M, 1);
       e_tilde = sqrt(sigma_e2)*white_noise;
       
    else
       % compute the Kronecker comb's period using the autocorrelation
       % method; ensure that the comb's frequency is higher than 200 Hz
       [r_e, lags] = xcorr(e);
       Tkc_min_samples = Tkc_min * Fs;
       
       r_e = r_e(lags >= 1);
       [pks, locs] = findpeaks(r_e);
       locs = locs(locs > Tkc_min_samples);
       Tkc = locs(1);
              
       % define the Kronecker comb
       Kcomb = zeros(1, M)';
       Kcomb(1:Tkc:end) = 1;

       
       % compute the excitation signal by scaling the Kronecker comb
       sigma_kc = (1/M)*sum(Kcomb.^2);
       sigma_g = sigma_e2 / sigma_kc;
       e_tilde = sigma_g .* Kcomb;
  
    end
    
    % synthetize the speech segment
    s_n = filter(1, a, e_tilde);   
    s_n = window .* s_n;
    
    % perform OLA
    out(i*Hop + 1:i*Hop + M) = out(i*Hop + 1 : i*Hop + M) + s_n;
    
    
end

% apply a 10 dB gain to the output signal (for playback volume adjustment)
gain_db = 10;
gain = 10^(gain_db/10);
out = gain*out;


% output audio playback
output_player = audioplayer(out, Fs);
play(output_player); 

% write the output file to disk 
filename = strcat(algorithm, "_output.wav");
audiowrite(filename, out, Fs);

% EOF