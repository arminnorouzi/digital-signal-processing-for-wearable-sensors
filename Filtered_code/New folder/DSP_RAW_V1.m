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
xlim([0 35])
subplot(3,2,2),plot(t,A_x_s2,'LineWidth', 1.5),legend('Sensor #2'),handlex9 = xlabel('time [s]');, handley9 = ylabel('a_x [m/s^2]');
xlim([0 35])
subplot(3,2,3),plot(t,A_y_s1,'LineWidth', 1.5),legend('Sensor #1'),handlex9 = xlabel('time [s]');, handley9 = ylabel('a_y [m/s^2]');
xlim([0 35])
subplot(3,2,4),plot(t,A_y_s2,'LineWidth', 1.5),legend('Sensor #2'),handlex9 = xlabel('time [s]');, handley9 = ylabel('a_y [m/s^2]');
xlim([0 35])
subplot(3,2,5),plot(t,A_z_s1,'LineWidth', 1.5),legend('Sensor #1'),handlex9 = xlabel('time [s]');, handley9 = ylabel('a_z [m/s^2]');
xlim([0 35])
subplot(3,2,6),plot(t,A_z_s2,'LineWidth', 1.5),legend('Sensor #2'),handlex9 = xlabel('time [s]');, handley9 = ylabel('a_z [m/s^2]');
xlim([0 35])
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

%% for slide
% h1= figure(30);
% set(gcf, 'Position', [0, 100, 1000, 550]);
% plot(f,P1A_x_s1f,f,P1A_x_s2f, f,P1A_y_s1f,f,P1A_y_s2f,f,P1A_z_s1f,f,P1A_z_s2f,'LineWidth', 1.5),legend('X direction of Sensor #1', 'X direction of Sensor #2', 'Y direction of Sensor #1', 'Y direction of Sensor #2', 'Z direction of Sensor #1', 'Z direction of Sensor #2'),handlex9 = xlabel('f (Hz)');, handley9 = ylabel('|P1(f)|');



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

n_ch = 5; % Chebyshev Type I order
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
subplot(2,3,1),plot(t, A_x_s1, t,A_x_s1_hnw, t,A_x_s1_hmw, t,A_x_s1_bw,'LineWidth', 1.5),title('X direction - sensor #1'),legend('Raw Data', 'Hanning Window','Hamming Window','Blackman Window'),handlex9 = xlabel('t [s]');, handley9 = ylabel('a_x [m/s^2]');
xlim([0 35])
subplot(2,3,2),plot(t, A_y_s1, t,A_y_s1_hnw, t,A_y_s1_hmw, t,A_y_s1_bw,'LineWidth', 1.5),title('y direction - sensor #1'),legend('Raw Data', 'Hanning Window','Hamming Window','Blackman Window'),handlex9 = xlabel('t [s]');, handley9 = ylabel('a_y [m/s^2]');
xlim([0 35])
subplot(2,3,3),plot(t, A_z_s1, t,A_z_s1_hnw, t,A_z_s1_hmw, t,A_z_s1_bw,'LineWidth', 1.5),title('z direction - sensor #1'),legend('Raw Data', 'Hanning Window','Hamming Window','Blackman Window'),handlex9 = xlabel('t [s]');, handley9 = ylabel('a_z [m/s^2]');
xlim([0 35])
subplot(2,3,4),plot(t, A_x_s1, t,A_x_s1_hnw, t,A_x_s1_hmw, t,A_x_s1_bw,'LineWidth', 1.5),title('X direction - sensor #1'),legend('Raw Data', 'Hanning Window','Hamming Window','Blackman Window'),handlex9 = xlabel('t [s]');, handley9 = ylabel('a_x [m/s^2]');
xlim([0 35])
subplot(2,3,5),plot(t, A_y_s1, t,A_y_s1_hnw, t,A_y_s1_hmw, t,A_y_s1_bw,'LineWidth', 1.5),title('y direction - sensor #1'),legend('Raw Data', 'Hanning Window','Hamming Window','Blackman Window'),handlex9 = xlabel('t [s]');, handley9 = ylabel('a_y [m/s^2]');
xlim([0 35])
subplot(2,3,6),plot(t, A_z_s1, t,A_z_s1_hnw, t,A_z_s1_hmw, t,A_z_s1_bw,'LineWidth', 1.5),title('z direction - sensor #1'),legend('Raw Data', 'Hanning Window','Hamming Window','Blackman Window'),handlex9 = xlabel('t [s]');, handley9 = ylabel('a_z [m/s^2]');
xlim([0 35])

figure(4)
set(gcf, 'Position', [0, 100, 1000, 550]);
subplot(2,3,1),plot(t, A_x_s2, t,A_x_s2_hnw, t,A_x_s2_hmw, t,A_x_s2_bw,'LineWidth', 1.5),title('X direction - sensor #2'),legend('Raw Data', 'Hanning Window','Hamming Window','Blackman Window'),handlex9 = xlabel('t [s]');, handley9 = ylabel('a_x [m/s^2]');
xlim([0 35])
subplot(2,3,2),plot(t, A_y_s2, t,A_y_s2_hnw, t,A_y_s2_hmw, t,A_y_s2_bw,'LineWidth', 1.5),title('y direction - sensor #2'),legend('Raw Data', 'Hanning Window','Hamming Window','Blackman Window'),handlex9 = xlabel('t [s]');, handley9 = ylabel('a_y [m/s^2]');
xlim([0 35])
subplot(2,3,3),plot(t, A_z_s2, t,A_z_s2_hnw, t,A_z_s2_hmw, t,A_z_s2_bw,'LineWidth', 1.5),title('z direction - sensor #2'),legend('Raw Data', 'Hanning Window','Hamming Window','Blackman Window'),handlex9 = xlabel('t [s]');, handley9 = ylabel('a_z [m/s^2]');
xlim([0 35])
subplot(2,3,4),plot(t, A_x_s2, t,A_x_s2_hnw, t,A_x_s2_hmw, t,A_x_s2_bw,'LineWidth', 1.5),title('X direction - sensor #2'),legend('Raw Data', 'Hanning Window','Hamming Window','Blackman Window'),handlex9 = xlabel('t [s]');, handley9 = ylabel('a_x [m/s^2]');
xlim([0 35])
subplot(2,3,5),plot(t, A_y_s2, t,A_y_s2_hnw, t,A_y_s2_hmw, t,A_y_s2_bw,'LineWidth', 1.5),title('y direction - sensor #2'),legend('Raw Data', 'Hanning Window','Hamming Window','Blackman Window'),handlex9 = xlabel('t [s]');, handley9 = ylabel('a_y [m/s^2]');
xlim([0 35])
subplot(2,3,6),plot(t, A_z_s2, t,A_z_s2_hnw, t,A_z_s2_hmw, t,A_z_s2_bw,'LineWidth', 1.5),title('z direction - sensor #2'),legend('Raw Data', 'Hanning Window','Hamming Window','Blackman Window'),handlex9 = xlabel('t [s]');, handley9 = ylabel('a_z [m/s^2]');
xlim([0 35])


%% FIR results For slides
% close all
% 
% figure(1)
% set(gcf, 'Position', [0, 100, 1000, 550]);
% subplot(2,1,1),plot(t, A_x_s1, t,A_x_s1_hnw, t,A_x_s1_hmw, t,A_x_s1_bw,'LineWidth', 1.5),title('X direction - sensor #1'),legend('Raw Data', 'Hanning Window','Hamming Window','Blackman Window'),handlex9 = xlabel('t [s]');, handley9 = ylabel('a_x [m/s^2]');
% xlim([10 30])
% subplot(2,1,2),plot(t, A_x_s1, t,A_x_s1_hnw, t,A_x_s1_hmw, t,A_x_s1_bw,'LineWidth', 1.5),title('X direction - sensor #1'),legend('Raw Data', 'Hanning Window','Hamming Window','Blackman Window'),handlex9 = xlabel('t [s]');, handley9 = ylabel('a_x [m/s^2]');
% xlim([10 30])
% 
% figure(2)
% set(gcf, 'Position', [0, 100, 1000, 550]);
% subplot(2,1,1),plot(t, A_y_s1, t,A_y_s1_hnw, t,A_y_s1_hmw, t,A_y_s1_bw,'LineWidth', 1.5),title('y direction - sensor #1'),legend('Raw Data', 'Hanning Window','Hamming Window','Blackman Window'),handlex9 = xlabel('t [s]');, handley9 = ylabel('a_y [m/s^2]');
% xlim([10 30])
% subplot(2,1,2),plot(t, A_y_s1, t,A_y_s1_hnw, t,A_y_s1_hmw, t,A_y_s1_bw,'LineWidth', 1.5),title('y direction - sensor #1'),legend('Raw Data', 'Hanning Window','Hamming Window','Blackman Window'),handlex9 = xlabel('t [s]');, handley9 = ylabel('a_y [m/s^2]');
% xlim([10 30])
% 
% figure(3)
% set(gcf, 'Position', [0, 100, 1000, 550]);
% subplot(2,1,1),plot(t, A_z_s1, t,A_z_s1_hnw, t,A_z_s1_hmw, t,A_z_s1_bw,'LineWidth', 1.5),title('z direction - sensor #1'),legend('Raw Data', 'Hanning Window','Hamming Window','Blackman Window'),handlex9 = xlabel('t [s]');, handley9 = ylabel('a_z [m/s^2]');
% xlim([10 30])
% subplot(2,1,2),plot(t, A_z_s1, t,A_z_s1_hnw, t,A_z_s1_hmw, t,A_z_s1_bw,'LineWidth', 1.5),title('z direction - sensor #1'),legend('Raw Data', 'Hanning Window','Hamming Window','Blackman Window'),handlex9 = xlabel('t [s]');, handley9 = ylabel('a_z [m/s^2]');
% xlim([10 30])
%% IIR Results            

figure(5)
set(gcf, 'Position', [0, 100, 1000, 550]);
subplot(2,3,1),plot(t, A_x_s1, t,A_x_s1_but, t,A_x_s1_ch,'LineWidth', 1.5),title('X direction - sensor #1'),legend('Raw Data','Butterworth','Chebyshev Type I'),handlex9 = xlabel('t [s]');, handley9 = ylabel('a_x [m/s^2]');
xlim([10 30])
subplot(2,3,2),plot(t, A_y_s1, t,A_y_s1_but, t,A_y_s1_ch,'LineWidth', 1.5),title('y direction - sensor #1'),legend('Raw Data','Butterworth','Chebyshev Type I'),handlex9 = xlabel('t [s]');, handley9 = ylabel('a_y [m/s^2]');
xlim([10 30])
subplot(2,3,3),plot(t, A_z_s1, t,A_z_s1_but, t,A_z_s1_ch,'LineWidth', 1.5),title('z direction - sensor #1'),legend('Raw Data','Butterworth','Chebyshev Type I'),handlex9 = xlabel('t [s]');, handley9 = ylabel('a_z [m/s^2]');
xlim([10 30])
subplot(2,3,4),plot(t, A_x_s1, t,A_x_s1_but, t,A_x_s1_ch,'LineWidth', 1.5),title('X direction - sensor #1'),legend('Raw Data','Butterworth','Chebyshev Type I'),handlex9 = xlabel('t [s]');, handley9 = ylabel('a_x [m/s^2]');
xlim([10 30])
subplot(2,3,5),plot(t, A_y_s1, t,A_y_s1_but, t,A_y_s1_ch,'LineWidth', 1.5),title('y direction - sensor #1'),legend('Raw Data','Butterworth','Chebyshev Type I'),handlex9 = xlabel('t [s]');, handley9 = ylabel('a_y [m/s^2]');
xlim([10 30])
subplot(2,3,6),plot(t, A_z_s1, t,A_z_s1_but, t,A_z_s1_ch,'LineWidth', 1.5),title('z direction - sensor #1'),legend('Raw Data','Butterworth','Chebyshev Type I'),handlex9 = xlabel('t [s]');, handley9 = ylabel('a_z [m/s^2]');
xlim([10 30])

figure(6)
set(gcf, 'Position', [0, 100, 1000, 550]);
subplot(2,3,1),plot(t, A_x_s2, t,A_x_s2_but, t,A_x_s2_ch,'LineWidth', 1.5),title('X direction - sensor #2'),legend('Raw Data','Butterworth','Chebyshev Type I'),handlex9 = xlabel('t [s]');, handley9 = ylabel('a_x [m/s^2]');
xlim([0 35])
subplot(2,3,2),plot(t, A_y_s2, t,A_y_s2_but, t,A_y_s2_ch,'LineWidth', 1.5),title('y direction - sensor #2'),legend('Raw Data','Butterworth','Chebyshev Type I'),handlex9 = xlabel('t [s]');, handley9 = ylabel('a_y [m/s^2]');
xlim([0 35])
subplot(2,3,3),plot(t, A_z_s2, t,A_z_s2_but, t,A_z_s2_ch,'LineWidth', 1.5),title('z direction - sensor #2'),legend('Raw Data','Butterworth','Chebyshev Type I'),handlex9 = xlabel('t [s]');, handley9 = ylabel('a_z [m/s^2]');
xlim([0 35])
subplot(2,3,4),plot(t, A_x_s2, t,A_x_s2_but, t,A_x_s2_ch,'LineWidth', 1.5),title('X direction - sensor #2'),legend('Raw Data','Butterworth','Chebyshev Type I'),handlex9 = xlabel('t [s]');, handley9 = ylabel('a_x [m/s^2]');
xlim([0 35])
subplot(2,3,5),plot(t, A_y_s2, t,A_y_s2_but, t,A_y_s2_ch,'LineWidth', 1.5),title('y direction - sensor #2'),legend('Raw Data','Butterworth','Chebyshev Type I'),handlex9 = xlabel('t [s]');, handley9 = ylabel('a_y [m/s^2]');
xlim([0 35])
subplot(2,3,6),plot(t, A_z_s2, t,A_z_s2_but, t,A_z_s2_ch,'LineWidth', 1.5),title('z direction - sensor #2'),legend('Raw Data','Butterworth','Chebyshev Type I'),handlex9 = xlabel('t [s]');, handley9 = ylabel('a_z [m/s^2]');
xlim([0 35])
%% IIR results For slides
% close all
% figure(1)
% set(gcf, 'Position', [0, 100, 1000, 550]);
% subplot(2,1,1),plot(t, A_x_s1, t,A_x_s1_but, t,A_x_s1_ch,'LineWidth', 1.5),title('X direction - sensor #1'),legend('Raw Data','Butterworth','Chebyshev Type I'),handlex9 = xlabel('t [s]');, handley9 = ylabel('a_x [m/s^2]');
% xlim([10 30])
% subplot(2,1,2),plot(t, A_x_s1, t,A_x_s1_but, t,A_x_s1_ch,'LineWidth', 1.5),title('X direction - sensor #1'),legend('Raw Data','Butterworth','Chebyshev Type I'),handlex9 = xlabel('t [s]');, handley9 = ylabel('a_x [m/s^2]');
% xlim([10 30])
% 
% figure(2)
% set(gcf, 'Position', [0, 100, 1000, 550]);
% subplot(2,1,1),plot(t, A_y_s1, t,A_y_s1_but, t,A_y_s1_ch,'LineWidth', 1.5),title('y direction - sensor #1'),legend('Raw Data','Butterworth','Chebyshev Type I'),handlex9 = xlabel('t [s]');, handley9 = ylabel('a_y [m/s^2]');
% xlim([10 30])
% subplot(2,1,2),plot(t, A_y_s1, t,A_y_s1_but, t,A_y_s1_ch,'LineWidth', 1.5),title('y direction - sensor #1'),legend('Raw Data','Butterworth','Chebyshev Type I'),handlex9 = xlabel('t [s]');, handley9 = ylabel('a_y [m/s^2]');
% xlim([10 30])
% 
% figure(3)
% set(gcf, 'Position', [0, 100, 1000, 550]);
% subplot(2,1,1),plot(t, A_z_s1, t,A_z_s1_but, t,A_z_s1_ch,'LineWidth', 1.5),title('z direction - sensor #1'),legend('Raw Data','Butterworth','Chebyshev Type I'),handlex9 = xlabel('t [s]');, handley9 = ylabel('a_z [m/s^2]');
% xlim([10 30])
% subplot(2,1,2),plot(t, A_z_s1, t,A_z_s1_but, t,A_z_s1_ch,'LineWidth', 1.5),title('z direction - sensor #1'),legend('Raw Data','Butterworth','Chebyshev Type I'),handlex9 = xlabel('t [s]');, handley9 = ylabel('a_z [m/s^2]');
% xlim([10 30])


%% Velocity Plot

%%%%%%%% Data input for integral %%%%%%%%
data1 = [A_x_s1 A_y_s1 A_z_s1]; %raw data sensor #1
data2 = [A_x_s2 A_y_s2 A_z_s2]; %raw data sensor #1
data1_but = [A_x_s1_but A_y_s1_but A_z_s1_but]; %Filtered Butterworth  sensor #1
data2_but = [A_x_s2_but A_y_s2_but A_z_s2_but]; %Filtered Butterworth  sensor #2

data1_ch = [A_x_s1_ch A_y_s1_ch A_z_s1_ch]; %Filtered Chebyshev Type I sensor #1
data2_ch = [A_x_s2_ch A_y_s2_ch A_z_s2_ch]; %Filtered Chebyshev Type I sensor #2

data1_hnw = [A_x_s1_hnw A_y_s1_hnw A_z_s1_hnw]; %Filtered hanning sensor #1
data2_hnw = [A_x_s2_hnw A_y_s2_hnw A_z_s2_hnw]; %Filtered hanning sensor #1

data1_hmw = [A_x_s1_hmw A_y_s1_hmw A_z_s1_hmw]; %Filtered hamming sensor #1
data2_hmw = [A_x_s2_hmw A_y_s2_hmw A_z_s2_hmw]; %Filtered hamming sensor #1

data1_bw = [A_x_s1_bw A_y_s1_bw A_z_s1_bw]; %Filtered Blackman sensor #1
data2_bw = [A_x_s2_bw A_y_s2_bw A_z_s2_bw]; %Filtered Blackman sensor #1

%%%%%%%% Velocity %%%%%%%%
%x
[time,vel_x]=velfinder(data1,data2,1);
[~,vel_x_but]=velfinder(data1_but,data2_but,1);
[~,vel_x_ch]=velfinder(data1_ch,data2_ch,1);
[~,vel_x_hnw]=velfinder(data1_hnw,data2_hnw,1);
[~,vel_x_hmw]=velfinder(data1_hmw,data2_hmw,1);
[~,vel_x_bw]=velfinder(data1_bw,data2_bw,1);

%y
[~,vel_y]=velfinder(data1,data2,2);
[~,vel_y_but]=velfinder(data1_but,data2_but,2);
[~,vel_y_ch]=velfinder(data1_ch,data2_ch,2);
[~,vel_y_hnw]=velfinder(data1_hnw,data2_hnw,2);
[~,vel_y_hmw]=velfinder(data1_hmw,data2_hmw,2);
[~,vel_y_bw]=velfinder(data1_bw,data2_bw,2);

%z
[~,vel_z]=velfinder(data1,data2,3);
[~,vel_z_but]=velfinder(data1_but,data2_but,3);
[~,vel_z_ch]=velfinder(data1_ch,data2_ch,3);
[~,vel_z_hnw]=velfinder(data1_hnw,data2_hnw,3);
[~,vel_z_hmw]=velfinder(data1_hmw,data2_hmw,3);
[~,vel_z_bw]=velfinder(data1_bw,data2_bw,3);


%%%%%%%% position %%%%%%%%
%x
[time_x,pos_x]=posfinder(data1,data2,1);
[~,pos_x_but]=posfinder(data1_but,data2_but,1);
[~,pos_x_ch]=posfinder(data1_ch,data2_ch,1);
[~,pos_x_hnw]=posfinder(data1_hnw,data2_hnw,1);
[~,pos_x_hmw]=posfinder(data1_hmw,data2_hmw,1);
[~,pos_x_bw]=posfinder(data1_bw,data2_bw,1);

%y
[~,pos_y]=posfinder(data1,data2,2);
[~,pos_y_but]=posfinder(data1_but,data2_but,2);
[~,pos_y_ch]=posfinder(data1_ch,data2_ch,2);
[~,pos_y_hnw]=posfinder(data1_hnw,data2_hnw,2);
[~,pos_y_hmw]=posfinder(data1_hmw,data2_hmw,2);
[~,pos_y_bw]=posfinder(data1_bw,data2_bw,2);

%z
[~,pos_z]=posfinder(data1,data2,3);
[~,pos_z_but]=posfinder(data1_but,data2_but,3);
[~,pos_z_ch]=posfinder(data1_ch,data2_ch,3);
[~,pos_z_hnw]=posfinder(data1_hnw,data2_hnw,3);
[~,pos_z_hmw]=posfinder(data1_hmw,data2_hmw,3);
[~,pos_z_bw]=posfinder(data1_bw,data2_bw,3);


%% Velocity

figure(7)
set(gcf, 'Position', [0, 100, 1000, 550]);
subplot(3,1,1),plot(time, vel_x(:,1), time, vel_x_hnw(:,1), time, vel_x_hmw(:,1), time, vel_x_bw(:,1),'LineWidth', 1.5),title('Velocity in X direction - sensor #1'),legend('Raw Data', 'Hanning Window','Hamming Window','Blackman Window'),handlex9 = xlabel('t [s]');, handley9 = ylabel('v_x [m/s]');
xlim([0 35])
subplot(3,1,2),plot(time, vel_y(:,1), time, vel_y_hnw(:,1), time, vel_y_hmw(:,1), time, vel_y_bw(:,1),'LineWidth', 1.5),title('Velocity in y direction - sensor #1'),legend('Raw Data', 'Hanning Window','Hamming Window','Blackman Window'),handlex9 = xlabel('t [s]');, handley9 = ylabel('v_y [m/s]');
xlim([0 35])
subplot(3,1,3),plot(time, vel_z(:,1), time, vel_z_hnw(:,1), time, vel_z_hmw(:,1), time, vel_z_bw(:,1),'LineWidth', 1.5),title('Velocity in z direction - sensor #1'),legend('Raw Data', 'Hanning Window','Hamming Window','Blackman Window'),handlex9 = xlabel('t [s]');, handley9 = ylabel('v_z [m/s]');
xlim([0 35])

figure(8)
set(gcf, 'Position', [0, 100, 1000, 550]);
subplot(3,1,1),plot(time, vel_x(:,2), time, vel_x_hnw(:,2), time, vel_x_hmw(:,2), time, vel_x_bw(:,2),'LineWidth', 1.5),title('Velocity in X direction - sensor #2'),legend('Raw Data', 'Hanning Window','Hamming Window','Blackman Window'),handlex9 = xlabel('t [s]');, handley9 = ylabel('v_x [m/s]');
xlim([0 35])
subplot(3,1,2),plot(time, vel_y(:,2), time, vel_y_hnw(:,2), time, vel_y_hmw(:,2), time, vel_y_bw(:,2),'LineWidth', 1.5),title('Velocity in y direction - sensor #2'),legend('Raw Data', 'Hanning Window','Hamming Window','Blackman Window'),handlex9 = xlabel('t [s]');, handley9 = ylabel('v_y [m/s]');
xlim([0 35])
subplot(3,1,3),plot(time, vel_z(:,2), time, vel_z_hnw(:,2), time, vel_z_hmw(:,2), time, vel_z_bw(:,2),'LineWidth', 1.5),title('Velocity in z direction - sensor #2'),legend('Raw Data', 'Hanning Window','Hamming Window','Blackman Window'),handlex9 = xlabel('t [s]');, handley9 = ylabel('v_z [m/s]');
xlim([0 35])

figure(9)
set(gcf, 'Position', [0, 100, 1000, 550]);
subplot(3,1,1),plot(time, vel_x(:,1), time, vel_x_but(:,1), time, vel_x_ch(:,1),'LineWidth', 1.5),title('Velocity in X direction - sensor #1'),legend('Raw Data','Butterworth','Chebyshev Type I'),handlex9 = xlabel('t [s]');, handley9 = ylabel('v_x [m/s]');
xlim([0 35])
subplot(3,1,2),plot(time, vel_y(:,1), time, vel_y_but(:,1), time, vel_y_ch(:,1),'LineWidth', 1.5),title('Velocity in y direction - sensor #1'),legend('Raw Data','Butterworth','Chebyshev Type I'),handlex9 = xlabel('t [s]');, handley9 = ylabel('v_y [m/s]');
xlim([0 35])
subplot(3,1,3),plot(time, vel_z(:,1), time, vel_z_but(:,1), time, vel_z_ch(:,1),'LineWidth', 1.5),title('Velocity in z direction - sensor #1'),legend('Raw Data','Butterworth','Chebyshev Type I'),handlex9 = xlabel('t [s]');, handley9 = ylabel('v_z [m/s]');
xlim([0 35])

figure(10)
set(gcf, 'Position', [0, 100, 1000, 550]);
subplot(3,1,1),plot(time, vel_x(:,2), time, vel_x_but(:,2), time, vel_x_ch(:,2),'LineWidth', 1.5),title('Velocity in X direction - sensor #2'),legend('Raw Data','Butterworth','Chebyshev Type I'),handlex9 = xlabel('t [s]');, handley9 = ylabel('v_x [m/s]');
xlim([0 35])
subplot(3,1,2),plot(time, vel_y(:,2), time, vel_y_but(:,2), time, vel_y_ch(:,2),'LineWidth', 1.5),title('Velocity in y direction - sensor #2'),legend('Raw Data','Butterworth','Chebyshev Type I'),handlex9 = xlabel('t [s]');, handley9 = ylabel('v_y [m/s]');
xlim([0 35])
subplot(3,1,3),plot(time, vel_z(:,2), time, vel_z_but(:,2), time, vel_z_ch(:,2),'LineWidth', 1.5),title('Velocity in z direction - sensor #2'),legend('Raw Data','Butterworth','Chebyshev Type I'),handlex9 = xlabel('t [s]');, handley9 = ylabel('v_z [m/s]');
xlim([0 35])

%%  Position results
figure(11)
set(gcf, 'Position', [0, 100, 1000, 550]);
subplot(3,1,1),plot(time_x, pos_x(:,1), time_x, pos_x_hnw(:,1), time_x, pos_x_hmw(:,1), time_x, pos_x_bw(:,1),'LineWidth', 1.5),title('Displacement in X direction - sensor #1'),legend('Raw Data', 'Hanning Window','Hamming Window','Blackman Window'),handlex9 = xlabel('t [s]');, handley9 = ylabel('x [m]');
xlim([0 35])
subplot(3,1,2),plot(time_x, pos_y(:,1), time_x, pos_y_hnw(:,1), time_x, pos_y_hmw(:,1), time_x, pos_y_bw(:,1),'LineWidth', 1.5),title('Displacement in y direction - sensor #1'),legend('Raw Data', 'Hanning Window','Hamming Window','Blackman Window'),handlex9 = xlabel('t [s]');, handley9 = ylabel('y [m]');
xlim([0 35])
subplot(3,1,3),plot(time_x, pos_z(:,1), time_x, pos_z_hnw(:,1), time_x, pos_z_hmw(:,1), time_x, pos_z_bw(:,1),'LineWidth', 1.5),title('Displacement in z direction - sensor #1'),legend('Raw Data', 'Hanning Window','Hamming Window','Blackman Window'),handlex9 = xlabel('t [s]');, handley9 = ylabel('z [m]');
xlim([0 35])

figure(12)
set(gcf, 'Position', [0, 100, 1000, 550]);
subplot(3,1,1),plot(time_x, pos_x(:,2), time_x, pos_x_hnw(:,2), time_x, pos_x_hmw(:,2), time_x, pos_x_bw(:,2),'LineWidth', 1.5),title('Displacement in X direction - sensor #2'),legend('Raw Data', 'Hanning Window','Hamming Window','Blackman Window'),handlex9 = xlabel('t [s]');, handley9 = ylabel('x [m]');
xlim([0 35])
subplot(3,1,2),plot(time_x, pos_y(:,2), time_x, pos_y_hnw(:,2), time_x, pos_y_hmw(:,2), time_x, pos_y_bw(:,2),'LineWidth', 1.5),title('Displacement in y direction - sensor #2'),legend('Raw Data', 'Hanning Window','Hamming Window','Blackman Window'),handlex9 = xlabel('t [s]');, handley9 = ylabel('y [m]');
xlim([0 35])
subplot(3,1,3),plot(time_x, pos_z(:,2), time_x, pos_z_hnw(:,2), time_x, pos_z_hmw(:,2), time_x, pos_z_bw(:,2),'LineWidth', 1.5),title('Displacement in z direction - sensor #2'),legend('Raw Data', 'Hanning Window','Hamming Window','Blackman Window'),handlex9 = xlabel('t [s]');, handley9 = ylabel('z [m]');


figure(13)
set(gcf, 'Position', [0, 100, 1000, 550]);
subplot(3,1,1),plot(time_x, pos_x(:,1), time_x, pos_x_but(:,1), time_x, pos_x_ch(:,1),'LineWidth', 1.5),title('Displacement in X direction - sensor #1'),legend('Raw Data','Butterworth','Chebyshev Type I'),handlex9 = xlabel('t [s]');, handley9 = ylabel('x [m]');
xlim([0 35])
subplot(3,1,2),plot(time_x, pos_y(:,1), time_x, pos_y_but(:,1), time_x, pos_y_ch(:,1),'LineWidth', 1.5),title('Displacement in y direction - sensor #1'),legend('Raw Data','Butterworth','Chebyshev Type I'),handlex9 = xlabel('t [s]');, handley9 = ylabel('y [m]');
xlim([0 35])
subplot(3,1,3),plot(time_x, pos_z(:,1), time_x, pos_z_but(:,1), time_x, pos_z_ch(:,1),'LineWidth', 1.5),title('Displacement in z direction - sensor #1'),legend('Raw Data','Butterworth','Chebyshev Type I'),handlex9 = xlabel('t [s]');, handley9 = ylabel('z [m]');
xlim([0 35])

figure(14)
set(gcf, 'Position', [0, 100, 1000, 550]);
subplot(3,1,1),plot(time_x, pos_x(:,2), time_x, pos_x_but(:,2), time_x, pos_x_ch(:,2),'LineWidth', 1.5),title('Displacement in X direction - sensor #2'),legend('Raw Data','Butterworth','Chebyshev Type I'),handlex9 = xlabel('t [s]');, handley9 = ylabel('x [m]');
xlim([0 35])
subplot(3,1,2),plot(time_x, pos_y(:,2), time_x, pos_y_but(:,2), time_x, pos_y_ch(:,2),'LineWidth', 1.5),title('Displacement in y direction - sensor #2'),legend('Raw Data','Butterworth','Chebyshev Type I'),handlex9 = xlabel('t [s]');, handley9 = ylabel('y [m]]');
xlim([0 35])
subplot(3,1,3),plot(time_x, pos_z(:,2), time_x, pos_z_but(:,2), time_x, pos_z_ch(:,2),'LineWidth', 1.5),title('Displacement in z direction - sensor #2'),legend('Raw Data','Butterworth','Chebyshev Type I'),handlex9 = xlabel('t [s]');, handley9 = ylabel('z [m]');
xlim([0 35])


%% Filter Config
figure(15)
set(gcf, 'Position', [0, 100, 1000, 550]);

% Pole-Zero map
subplot(3,5,1),zplane(b_hann,a_hann),title('Hanning Window');
subplot(3,5,2),zplane(b_hamm,a_hamm),title('Hamming Window');
subplot(3,5,3),zplane(b_bl,a_bl),title('Blackman Window');
subplot(3,5,4),zplane(b_bw,a_bw),title('Butterworth');
subplot(3,5,5),zplane(b_ch,a_ch),title('Chebyshev Type I');

% Impulse Response
subplot(3,5,6),impz(b_hann,a_hann),title('Hanning Window');
subplot(3,5,7),impz(b_hamm,a_hamm),title('Hamming Window');
subplot(3,5,8),impz(b_bl,a_bl),title('Blackman Window');
subplot(3,5,9),impz(b_bw,a_bw),title('Butterworth');
subplot(3,5,10),impz(b_ch,a_ch),title('Chebyshev Type I');

% Step Response
subplot(3,5,11),stepz(b_hann,a_hann),title('Hanning Window');
subplot(3,5,12),stepz(b_hamm,a_hamm),title('Hamming Window');
subplot(3,5,13),stepz(b_bl,a_bl),title('Blackman Window');
subplot(3,5,14),stepz(b_bw,a_bw),title('Butterworth');
subplot(3,5,15),stepz(b_ch,a_ch),title('Chebyshev Type I');

%% Filter Shape

% Hanning Window
figure(16)
set(gcf, 'Position', [10, 200, 1200, 800]);
freqz(b_hann,a_hann),title('Hanning Window');

% Hamming Window
figure(17)
set(gcf, 'Position', [10, 200, 1200, 800]);
freqz(b_hamm,a_hamm),title('Hamming Window');

% Blackman Window
figure(18)
set(gcf, 'Position', [10, 200, 1200, 800]);
freqz(b_bl,a_bl),title('Blackman Window');

% Butterworth
figure(19)
set(gcf, 'Position', [10, 200, 1200, 800]);
subplot(2,5,[4,9]),freqz(b_bw,a_bw),title('Butterworth');

% Chebyshev Type I
figure(20)
set(gcf, 'Position', [10, 200, 1200, 800]);
subplot(2,5,[5,10]),freqz(b_ch,a_ch),title('Chebyshev Type I');


%% For slide
figure(21)
findpeaks(vel_z(:,1),'MinPeakDistance',100);title('Velocity in z direction - sensor #2'),legend('Raw Data','Butterworth','Chebyshev Type I'),handlex9 = xlabel('t [s]');, handley9 = ylabel('v_z [m/s]');
xlim([1400 2400])