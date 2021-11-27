%----------------------------------------------%
%           *** DAAP - HOMEWORK #2 ***         %
%----------------------------------------------%
%                DOA estimation                %
%----------------------------------------------%
% Giovanni Affatato, Roberto Alessandri        %
%----------------------------------------------%

clear; close all; clc;

%% Parameters

% import the array vector y
load("array_data_64_mics.mat");
%load("array_data_8_mics.mat")

% sampling frequency
Fs = 8000; 
% speed of sound [m/s]
c = 340;    
% number of sources
N_src = 2;

% M: number of microphones
% N: microphone signal length
[M, N] = size(y);

% determine the distance between two mics in meters (see anti-aliasing condition)
FN = Fs/2; % Nyquist
lambda_min = c/FN;
d = lambda_min/2;

%% Frequency estimation

% plot the magnitude spectrum of a microphone signal
y_1 = y(1, :);
y_1_Magnitude = abs(y_1);
freq_ax = 0:Fs/N:Fs*(1-1/N);

figure(1)
plot(freq_ax, y_1_Magnitude)
xlabel('Frequency [Hz]')
ylabel('Magnitude')
title('Magnitude spectrum of a microphone signal')
saveas(gcf,'Figure1.fig')

% estimate of the angular frequency of source signal 1
y1_MagNyq = y_1_Magnitude(1:FN*N/Fs);
[peaks, locs] = findpeaks(y1_MagNyq);
[peaks_sorted, indexes] = sort(peaks,'descend');
f_1 = (locs(indexes(1))-1)*Fs/N;
w_c1 = f_1*2*pi;


% estimate of the angular frequency of source signal 2
f_2 = (locs(indexes(2))-1)*Fs/N;
w_c2 = f_2*2*pi;

% the array of the estimated angular frequencies
w_c = [w_c1, w_c2];

%% Delay-and-sum beamformer

% the angles to be examined (candidates for our DOA estimate)
angles = (-pi/2:pi/180:pi/2); % in radiants

% the number of the angles to be examined
N_angles = length(angles);

% pre-allocate the DAS pseudo-spectrum
DAS_pseudo = zeros(N_angles,1);
% pre-allocate the DAS DOA estimates of both sources 
DAS_DOA = zeros(N_src,1);
% pre-allocate the steering vectors for both sources
a = zeros(M, N_angles, N_src);
% pre-allocate the spatial filter coefficients for both sources
h_DOA = zeros(M, N_src);

% sample estimate of the covariance matrix of the array vector
R = (1/(N-1)) * y * (y');

for i=1:N_src % for each source...
    
    for j=1:N_angles % for every angle to be examined...
        
        % compute the spatial frequency
        w_s = w_c(i) * d * sin(angles(j)) / c;%
        % compute the steering vector for the j-th angle and i-th source 
        a(:,j,i) = exp(-1i*w_s.*(0:M-1));%
        % save the value of the pseudo-spectrum for the j-th angle
        DAS_pseudo(j) = (a(:,j,i)' * R * a(:,j,i)) / (M^2);
        
    end
    
    % remove the spurious imaginary part
    DAS_pseudo = real(DAS_pseudo);
    
    % plot the DAS pseudo-spectrum
    figure(i+1)
    subplot(211)
    plot(rad2deg(angles), DAS_pseudo)
    xlabel('Angle [deg]')
    ylabel('Pseudo-spectrum')
    title("Delay-and-sum beamformer: pseudo-spectrum at " + num2str(w_c(i)/(2*pi)) + " Hz")
    
    % estimate the DOA of the i-th source as the angle for which the 
    % DAS pseudo-spectrum exhibits the most prominent peak 
    [~, DAS_DOA(i)] = max(DAS_pseudo);
    DAS_DOA(i) = angles(DAS_DOA(i));
    
    % compute the spatial frequency associated to the estimated DOA
    w_s_DOA = w_c(i)*d*sin(DAS_DOA(i))/c;
    % compute the steering vector associated to the estimated DOA
    a_DOA = exp(-1i*w_s_DOA.*(0:M-1));%
    % compute and store the spatial filter associated to the estimated DOA
    h_DOA(:,i) = a_DOA/M;
    
    % check for h_DOA denominator = a_DOA*a_DOA' = M
    M__ = a_DOA * a_DOA';
    
     
    % compute the spatial response using the spatial filter associated to
    % the estimated DOA and the pre-computed steering vectors 
    h_DOA_H_i = h_DOA(:,i)';
    DAS_spatial_response = h_DOA_H_i*a(:, :, i);%
    
    % remove the spurious imaginary part
    DAS_spatial_response = real(DAS_spatial_response);
    
    % plot the DAS spatial response
    figure(i+1)
    subplot(212)
    polarplot(angles, DAS_spatial_response);
    title("Delay-and-sum beamformer: beam pattern at " + num2str(w_c(i)/(2*pi)) + " Hz")
    thetalim([-90, 90])
    saveas(gcf,strcat('Figure', int2str(i+1)));
    
end

% Display the DAS DOA estimates (in degrees)
disp('DELAY-AND-SUM:')
disp(rad2deg(DAS_DOA(1)))
disp(rad2deg(DAS_DOA(2)))

%% Apply spatial filtering to the array vector
s_hat_1 = h_DOA(:,1)'*y;%
s_hat_2 = h_DOA(:,2)'*y;%

%% Play a microphone signal and the two beamformer outputs
pause
%% play a microphone signal (contains a mixture of the two sources)
soundsc(real(ifft(y_1)), Fs)
pause
%% play the filtered singal (beamformer steered towards source 1)
soundsc(real(ifft(s_hat_1)), Fs)
pause
%% play the filtered singal (beamformer steered towards source 2)
soundsc(real(ifft(s_hat_2)), Fs)

pause

%% Plot the pseudo-spectrum as a function of both frequency and angle of arrival

% the angular frequencies to be examined
angular_freqs = (0:10:FN)*2*pi;
% the number of angular frequencies to be examined
N_freqs = length(angular_freqs);

% pre-allocate the DAS pseudo-spectrum
pseudo_ = zeros(N_freqs, N_angles);

for k=1:N_freqs % for every angular frequency to be examined...
    for j=1:N_angles % for every angle to be examined...
        % compute the spatial frequency
        w_s_ = angular_freqs(k)*d*sin(angles(j))/c;%
        % compute the steering vector
        a_ = conj(exp(-1i*w_s_.*(0:M-1)))';%
        % compute and store the value of the pseudo-spectrum for the given frequency and angle
        pseudo_(k,j) = (a_' * R * a_) / (M^2);
    end
end

% remove the spurious imaginary part
pseudo_ = real(pseudo_);

% plot the DAS pseudo-spectrum as a function of frequency (in Hz) and the 
% angles of arrival (in degrees)
figure(4)
imagesc(rad2deg(angles), angular_freqs/(2*pi), pseudo_)
xlabel('Angle [deg]')
ylabel('Frequency [Hz]')
set(gca, 'YDir', 'normal')
colorbar()
title('Delay-and-sum beamformer: pseudo-spectrum')
saveas(gcf,'Figure4.fig')

%% Parametric methods

% eigenvalue decomposition of the covariance matrix R
[Q, L] = eig(R);
e2 = norm(R*Q-Q*L);

% sort the eigenvalues in descending order
[l,ind] = sort(diag(L), 'descend');
L = L(ind,ind);

% permute the columns of the eigenvector matrix accordingly
Q = Q(:,ind);

% e = 0 if Q sorted and L sorted are again a eigendecomposition of R
e1 = norm(R*Q-Q*L);
e = abs(e1 - e2); %OK

%% MUSIC

% compute the matrix whose columns span the noise subspace
V = Q(:, N_src + 1:M);

% pre-allocate the MUSIC pseudo-spectrum
MUSIC_pseudo = zeros(N_angles,1);

% pre-allocate the DOA estimates of both sources
MUSIC_DOA = zeros(N_src,1);

for i = 1:N_src % for each source...
    
    for j = 1:N_angles 
        % compute the MUSIC pseudo-spectrum using the pre-computed steering vectors
        MUSIC_pseudo(j) = 1 / (a(:,j,i)' * (V * V') * a(:,j,i));
    end
    
    % remove the spurious imaginary part
    MUSIC_pseudo = real(MUSIC_pseudo);
        
    % estimate the DOA of the i-th source as the angle for which the 
    % MUSIC pseudo-spectrum exhibits the most prominent peak
    [~, MUSIC_DOA(i)] = max(MUSIC_pseudo);
    MUSIC_DOA(i) = angles(MUSIC_DOA(i));
        
    % plot the MUSIC pseudo-spectrum
    figure(5)
    subplot(2,1,i)
    plot(rad2deg(angles), MUSIC_pseudo)
    xlabel('Angle [deg]')
    ylabel('Pseudo-spectrum')
    title("MUSIC pseudo-spectrum at " + num2str(w_c(i)/(2*pi)) + " Hz")
    
end
saveas(gcf,'Figure5.fig')

% display the MUSIC DOA estimates (in degrees)
disp('MUSIC:')
disp(rad2deg(MUSIC_DOA(1)))
disp(rad2deg(MUSIC_DOA(2)))

%% ESPRIT

% compute the matrix whose columns span the source subspace
U = Q(:, 1: N_src);

% compute the matrices of eigenvectors associated to the sub-arrays
M_ = M - 1;
U1 = [eye(M_) zeros(M_, 1)] * U;
U2 = [zeros(M_, 1) eye(M_)] * U;

% compute the least-square estimate of the matrix Phi
Phi_hat = pinv(U1) * U2;

% compute the eigenvalues of the matrix Phi_hat
nu = eig(Phi_hat);

% estimate the spatial frequencies
w_sl(1) = -angle(nu(1));
w_sl(2) = -angle(nu(2));

% estimate the DOAs
ESPRIT_DOA_1 = asin((w_sl(1) * c) / (w_c(1) * d));
ESPRIT_DOA_2 = asin((w_sl(2) * c) / (w_c(2) * d));

% display the ESPRIT DOA estimates (in degrees)
disp('ESPRIT:')
disp(rad2deg(ESPRIT_DOA_1))
disp(rad2deg(ESPRIT_DOA_2))