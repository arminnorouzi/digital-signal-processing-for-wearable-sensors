clc
clear all
close all

%% Import Data
s1 = importdata('S1.mat');
s2 = importdata('S2.mat');

%% Data classification

ctr = input('Acceleration  type (1= Acceleration, 2 = Free acceleration) =      ');

    if ctr ==1 
        n1 = 3; n2 = 4; n3 = 5;
    else
        n1 = 6; n2 = 7; n3 = 8;
    end
A_x_s1 = s1(:,n1);
A_y_s1 = s1(:,n2);
A_z_s1 = s1(:,n3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A_x_s2 = s2(:,n1);
A_y_s2 = s2(:,n2);
A_z_s2 = s2(:,n3);

%% Finding Time
Fs = 100;            % Sampling frequency                    
T = 1/Fs;             % Sampling period       
L = size(A_z_s1,1);             % Length of signal
t = (0:L-1)*T;        % Time vector
f = Fs*(0:(L/2))/L;

%% plot - time domain
h1= figure(1),set(gcf, 'Position', [0, 100, 1000, 550]);
subplot(3,2,1),plot(t,A_x_s1,'LineWidth', 1.5),legend('Sensor #1'),handlex9 = xlabel('time [s]');, handley9 = ylabel('a_x [m/s^2]');
subplot(3,2,2),plot(t,A_x_s2,'LineWidth', 1.5),legend('Sensor #2'),handlex9 = xlabel('time [s]');, handley9 = ylabel('a_x [m/s^2]');
subplot(3,2,3),plot(t,A_y_s1,'LineWidth', 1.5),legend('Sensor #1'),handlex9 = xlabel('time [s]');, handley9 = ylabel('a_y [m/s^2]');
subplot(3,2,4),plot(t,A_y_s2,'LineWidth', 1.5),legend('Sensor #2'),handlex9 = xlabel('time [s]');, handley9 = ylabel('a_y [m/s^2]');
subplot(3,2,5),plot(t,A_z_s1,'LineWidth', 1.5),legend('Sensor #1'),handlex9 = xlabel('time [s]');, handley9 = ylabel('a_z [m/s^2]');
subplot(3,2,6),plot(t,A_z_s2,'LineWidth', 1.5),legend('Sensor #2'),handlex9 = xlabel('time [s]');, handley9 = ylabel('a_z [m/s^2]');

%% Freq Domain
A_x_s1f = fft(A_x_s1);
P2A_x_s1f = abs(A_x_s1f/L);
P1A_x_s1f = P2A_x_s1f(1:L/2+1);
P1A_x_s1f(2:end-1) = 2*P1A_x_s1f(2:end-1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A_y_s1f = fft(A_y_s1);
P2A_y_s1f = abs(A_y_s1f/L);
P1A_y_s1f = P2A_y_s1f(1:L/2+1);
P1A_y_s1f(2:end-1) = 2*P1A_y_s1f(2:end-1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A_z_s1f = fft(abs(A_z_s1));
P2A_z_s1f = abs(A_z_s1f/L);
P1A_z_s1f = P2A_z_s1f(1:L/2+1);
P1A_z_s1f(2:end-1) = 2*P1A_z_s1f(2:end-1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A_x_s2f = fft(A_x_s2);
P2A_x_s2f = abs(A_x_s2f/L);
P1A_x_s2f = P2A_x_s2f(1:L/2+1);
P1A_x_s2f(2:end-1) = 2*P1A_x_s2f(2:end-1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A_y_s2f = fft(A_y_s2);
P2A_y_s2f = abs(A_y_s2f/L);
P1A_y_s2f = P2A_y_s2f(1:L/2+1);
P1A_y_s2f(2:end-1) = 2*P1A_y_s2f(2:end-1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A_z_s2f = fft(A_z_s2);
P2A_z_s2f = abs(A_z_s2f/L);
P1A_z_s2f = P2A_z_s2f(1:L/2+1);
P1A_z_s2f(2:end-1) = 2*P1A_z_s2f(2:end-1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h1= figure(2);
set(gcf, 'Position', [0, 100, 1000, 550]);
subplot(3,1,1),plot(f,P1A_x_s1f,f,P1A_x_s2f,'LineWidth', 1.5),title('X direction'),legend('Sensor #1', 'Sensor #2'),handlex9 = xlabel('f (Hz)');, handley9 = ylabel('|P1(f)|');
subplot(3,1,2),plot(f,P1A_y_s1f,f,P1A_y_s2f,'LineWidth', 1.5),title('Y direction'),legend('Sensor #1', 'Sensor #2'),handlex9 = xlabel('f (Hz)');, handley9 = ylabel('|P1(f)|');
subplot(3,1,3),plot(f,P1A_z_s1f,f,P1A_z_s2f,'LineWidth', 1.5),title('Z direction'),legend('Sensor #1', 'Sensor #2'),handlex9 = xlabel('f (Hz)');, handley9 = ylabel('|P1(f)|');


%% Filter Design
fs = 100;

%%%%% FIR %%%%
PB = 10; % Pass Band Edge (Hz)
SB = 20; % Stop Band Edge (Hz)

TW = (SB-PB);
f1 = PB+(TW/2);
omega = 2*pi*(f1/fs);


%%%%%%%% Hanning Window %%%%%%%%

Nhann = 3.32*(fs/TW); % Hanning Window
Nhann = round(Nhann);

if rem(Nhann,2)==0
 Nhann = Nhann+1;
end

n_hann = [-(Nhann-1)/2:(Nhann-1)/2];
h_hann1 = (omega/pi).*sinc((omega/pi)*n_hann);

w_hann = hann(Nhann)';
h_hann = h_hann1.*w_hann;

b_hann = h_hann;
a_hann = 1;



A_x_s1_hnw = filtfilt(b_hann,a_hann,A_x_s1);
A_y_s1_hnw = filtfilt(b_hann,a_hann,A_y_s1);
A_z_s1_hnw = filtfilt(b_hann,a_hann,A_z_s1);

A_x_s2_hnw = filtfilt(b_hann,a_hann,A_x_s2);
A_y_s2_hnw = filtfilt(b_hann,a_hann,A_y_s2);
A_z_s2_hnw = filtfilt(b_hann,a_hann,A_z_s2);




%%%%%%%% Hamming Window %%%%%%%%
Nhamm = 3.44*(fs/TW); % Hamming Window
Nhamm = round(Nhamm);

if rem(Nhamm,2)==0
 Nhamm = Nhamm+1;
end

n_hamm = [-(Nhamm-1)/2:(Nhamm-1)/2];
h_hamm1 = (omega/pi).*sinc((omega/pi)*n_hamm);

w_hamm = hamming(Nhamm)';
h_hamm = h_hamm1.*w_hamm;

b_hamm = h_hamm;
a_hamm = 1;

A_x_s1_hmw = filtfilt(b_hamm,a_hamm,A_x_s1);
A_y_s1_hmw = filtfilt(b_hamm,a_hamm,A_y_s1);
A_z_s1_hmw = filtfilt(b_hamm,a_hamm,A_z_s1);

A_x_s2_hmw = filtfilt(b_hamm,a_hamm,A_x_s2);
A_y_s2_hmw = filtfilt(b_hamm,a_hamm,A_y_s2);
A_z_s2_hmw = filtfilt(b_hamm,a_hamm,A_z_s2);


%%%%%%%% Blackman Window %%%%%%%%
Nbl = 5.98*(fs/TW); % Blackman Window
Nbl = round(Nbl);

if rem(Nbl,2)==0
 Nbl = Nbl+1;
end

n_bl = [-(Nbl-1)/2:(Nbl-1)/2];
h_bl1 = (omega/pi).*sinc((omega/pi)*n_bl);

w_bl = blackman(Nbl)';
h_bl = h_bl1.*w_bl;

b_bl = h_bl;
a_bl = 1;

A_x_s1_bw = filtfilt(b_bl,a_bl,A_x_s1);
A_y_s1_bw = filtfilt(b_bl,a_bl,A_y_s1);
A_z_s1_bw = filtfilt(b_bl,a_bl,A_z_s1);

A_x_s2_bw = filtfilt(b_bl,a_bl,A_x_s2);
A_y_s2_bw = filtfilt(b_bl,a_bl,A_y_s2);
A_z_s2_bw = filtfilt(b_bl,a_bl,A_z_s2);

%%% IIR Filter %%%%

%%%%%%%% Butterworth %%%%%%%%
fc = 10; % cut-off freq
fs = 100; % sampling freq
n_bw = 4; % order of filter
[b_bw,a_bw] = butter(n_bw,fc/(fs/2));

A_x_s1_but = filtfilt(b_bw,a_bw,A_x_s1);
A_y_s1_but = filtfilt(b_bw,a_bw,A_y_s1);
A_z_s1_but = filtfilt(b_bw,a_bw,A_z_s1);

A_x_s2_but = filtfilt(b_bw,a_bw,A_x_s2);
A_y_s2_but = filtfilt(b_bw,a_bw,A_y_s2);
A_z_s2_but = filtfilt(b_bw,a_bw,A_z_s2);


%%%%%%%% Chebyshev Type I %%%%%%%%

n_ch = 4; % Chebyshev Type I order
fp_ch = 10; % Pass band edge frequency (Hz)
Rp_ch = 0.001; % Pass band ripple (dB)
wp_ch = (2*pi*fp_ch)/(fs*pi);
[b_ch,a_ch] = cheby1(n_ch,Rp_ch,wp_ch);  

A_x_s1_ch = filtfilt(b_ch,a_ch,A_x_s1);
A_y_s1_ch = filtfilt(b_ch,a_ch,A_y_s1);
A_z_s1_ch = filtfilt(b_ch,a_ch,A_z_s1);

A_x_s2_ch = filtfilt(b_ch,a_ch,A_x_s2);
A_y_s2_ch = filtfilt(b_ch,a_ch,A_y_s2);
A_z_s2_ch = filtfilt(b_ch,a_ch,A_z_s2);


%%%%%%%% Raw Data vs Filtered Data %%%%%%%%

%% FIR results
figure(3)
set(gcf, 'Position', [0, 100, 1000, 550]);
subplot(3,1,1),plot(t, A_x_s1, t,A_x_s1_hnw, t,A_x_s1_hmw, t,A_x_s1_bw,'LineWidth', 1.5),title('X direction - sensor #1'),legend('Raw Data', 'Hanning Window','Hamming Window','Blackman Window'),handlex9 = xlabel('t [s]');, handley9 = ylabel('a_x [m/s^2]');
subplot(3,1,2),plot(t, A_y_s1, t,A_y_s1_hnw, t,A_y_s1_hmw, t,A_y_s1_bw,'LineWidth', 1.5),title('y direction - sensor #1'),legend('Raw Data', 'Hanning Window','Hamming Window','Blackman Window'),handlex9 = xlabel('t [s]');, handley9 = ylabel('a_x [m/s^2]');
subplot(3,1,3),plot(t, A_z_s1, t,A_z_s1_hnw, t,A_z_s1_hmw, t,A_z_s1_bw,'LineWidth', 1.5),title('z direction - sensor #1'),legend('Raw Data', 'Hanning Window','Hamming Window','Blackman Window'),handlex9 = xlabel('t [s]');, handley9 = ylabel('a_x [m/s^2]');

figure(4)
set(gcf, 'Position', [0, 100, 1000, 550]);
subplot(3,1,1),plot(t, A_x_s2, t,A_x_s2_hnw, t,A_x_s2_hmw, t,A_x_s2_bw,'LineWidth', 1.5),title('X direction - sensor #2'),legend('Raw Data', 'Hanning Window','Hamming Window','Blackman Window'),handlex9 = xlabel('t [s]');, handley9 = ylabel('a_x [m/s^2]');
subplot(3,1,2),plot(t, A_y_s2, t,A_y_s2_hnw, t,A_y_s2_hmw, t,A_y_s2_bw,'LineWidth', 1.5),title('y direction - sensor #2'),legend('Raw Data', 'Hanning Window','Hamming Window','Blackman Window'),handlex9 = xlabel('t [s]');, handley9 = ylabel('a_x [m/s^2]');
subplot(3,1,3),plot(t, A_z_s2, t,A_z_s2_hnw, t,A_z_s2_hmw, t,A_z_s2_bw,'LineWidth', 1.5),title('z direction - sensor #2'),legend('Raw Data', 'Hanning Window','Hamming Window','Blackman Window'),handlex9 = xlabel('t [s]');, handley9 = ylabel('a_x [m/s^2]');


%% IIR Results

figure(5)
set(gcf, 'Position', [0, 100, 1000, 550]);
subplot(3,1,1),plot(t, A_x_s1, t,A_x_s1_but, t,A_x_s1_ch,'LineWidth', 1.5),title('X direction - sensor #1'),legend('Raw Data','Butterworth','Chebyshev Type I'),handlex9 = xlabel('t [s]');, handley9 = ylabel('a_x [m/s^2]');
subplot(3,1,2),plot(t, A_y_s1, t,A_y_s1_but, t,A_y_s1_ch,'LineWidth', 1.5),title('y direction - sensor #1'),legend('Raw Data','Butterworth','Chebyshev Type I'),handlex9 = xlabel('t [s]');, handley9 = ylabel('a_x [m/s^2]');
subplot(3,1,3),plot(t, A_z_s1, t,A_z_s1_but, t,A_z_s1_ch,'LineWidth', 1.5),title('z direction - sensor #1'),legend('Raw Data','Butterworth','Chebyshev Type I'),handlex9 = xlabel('t [s]');, handley9 = ylabel('a_x [m/s^2]');

figure(6)
set(gcf, 'Position', [0, 100, 1000, 550]);
subplot(3,1,1),plot(t, A_x_s2, t,A_x_s2_but, t,A_x_s2_ch,'LineWidth', 1.5),title('X direction - sensor #2'),legend('Raw Data','Butterworth','Chebyshev Type I'),handlex9 = xlabel('t [s]');, handley9 = ylabel('a_x [m/s^2]');
subplot(3,1,2),plot(t, A_y_s2, t,A_y_s2_but, t,A_y_s2_ch,'LineWidth', 1.5),title('y direction - sensor #2'),legend('Raw Data','Butterworth','Chebyshev Type I'),handlex9 = xlabel('t [s]');, handley9 = ylabel('a_x [m/s^2]');
subplot(3,1,3),plot(t, A_z_s2, t,A_z_s2_but, t,A_z_s2_ch,'LineWidth', 1.5),title('z direction - sensor #2'),legend('Raw Data','Butterworth','Chebyshev Type I'),handlex9 = xlabel('t [s]');, handley9 = ylabel('a_x [m/s^2]');




