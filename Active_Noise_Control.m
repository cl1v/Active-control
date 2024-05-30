%% Setup Workspace
clear; clc
load('Plant_Response.mat')
load('Primary_Response')

%%  <<<<<<<<<<<<<<       PART 1  Plant Response Identification    >>>>>>>>>

% Calculate the plant response, between the loudspeaker and error
% microphone using the FRF estimator

nfft = 2^10;        % set number of points in FFT
noverlap = nfft/2;  % set overlap
win = nfft;         % set length of default hamming window

[Sxx, ~]=pwelch(u_ID, win, noverlap, nfft, Fs);     % power spectral density (PSD) of input signal
[Sxy, f]=cpsd(e_ID, u_ID, win, noverlap, nfft, Fs); % cross PSD between signal mneasured at error sensor and identification signal driving the control loudspeaker
G = Sxy./Sxx;                                       % plant transfer response calculation  - H1 estimator

[G,f] = tfestimate(u_ID,e_ID,win,noverlap,f,Fs);

% Plot the magnitude and phase response of the plant
figure()
subplot(2,1,1)
semilogx(f, 20*log10(abs(G)))
xlim([0, max(f)])
xlabel('Frequency, Hz')
ylabel('20 log_{10}|G|, dB')
subplot(2,1,2)
semilogx(f, (angle(G)))
xlim([0, max(f)])
xlabel('Frequency, Hz')
ylabel('\angle G, rad')

%% Calculate the impulse response
I = 128;                                         % Plant model length
G_ifft=fftshift(ifft([G; conj(G(end:-1:2))]));  % inverse FFT frequency domain plant response - operating on conjungate symmetric response between 0 and Fs
[G_FIR,~]=invfreqz(G,f/max(f)*pi,I-1,0);        % alternative means to calculate impulse response
t = 0:1/Fs:(I-1)/Fs;
figure()
plot(t, G_FIR,'k', t, G_ifft((length(G_ifft)-1)/2+1:1:(length(G_ifft)-1)/2+I),'r--')
xlim([0, max(t)])
xlabel('t, sec')
ylabel('G_{FIR}')
%%

G_model=freqz(G_FIR,1,f,Fs);
figure()
subplot(2,1,1)
semilogx(f, 20*log10(abs(G)),'k', f, 20*log10(abs(G_model)),'r')
subplot(2,1,2)
semilogx(f, (angle(G)),'k', f, (angle(G_model)),'r')

%%  <<<<<<<<<<<<<<       PART 2    Optimal Harmonic Controller  >>>>>>>>>>>
beta = 1e-2;                               % Regularisation
% Calculate one-sided fft of the disturbance signal

d_k = fft(d, 2*(length(f)-1));              % FFT of disturbance time history
d_k = (d_k(1:length(f))/(2*length(f))).';   % truncate to give the one-sided spectrum
u_k = -d_k./(G+beta);                       % caluclate the optimal control signal at each frequency
e_k = d_k + G.*u_k;                         % calculate the error signal at each frequency

% Plot
figure()
semilogx(f, 20*log10(abs(d_k)),'k', f, 20*log10(abs(e_k)), 'r')
xlim([0, max(f)])
xlabel('Frequency, Hz')
ylabel('20log_{10}(E_k), dB')


%%  <<<<<<<<<<<<       PART 3  Unconstrianed Broadband Control    >>>>>>>>>

S_xd = cpsd(x, d, win, noverlap, nfft, Fs);     % CPSD between the reference and disturbance signals
S_xx = pwelch(x, win, noverlap, nfft, Fs);      % PSD of the reference signal
S_dd = pwelch(d, win, noverlap, nfft, Fs);      % PSD of the disturbance signal
W_opt = -S_xd./(G.*S_xx);                       % calculate the optimal, unconstrained broadband controller
S_ee = S_dd-(abs(S_xd).^2)./S_xx;               % power spectrum of error signal (as in handout)

% Plot Optimal Filter Response
figure()
subplot(2,1,1)
plot(f, 20*log10(abs(W_opt)))
xlabel('Frequency, Hz')
ylabel('20 log_{10}|W_{opt}| dB')
subplot(2,1,2)
plot(f, (angle(W_opt)))
xlabel('Frequency, Hz')
ylabel('\angle W_{opt}, rad')

%% Plot control performance
figure()
semilogx(f, 10*log10(S_dd),'k', f, 10*log10(S_ee), 'r')
xlabel('Frequency, Hz')
xlim([0, max(f)])
ylabel('20log_{10}(S_{ee}), dB')
legend('Disturbance', 'Error')
%% Coherence

C_xd = mscohere(x, d, win, noverlap, nfft, Fs);     % CPSD between the reference and disturbance signals
%%
figure()
semilogx(f, C_xd)
xlabel('Frequency, Hz')
ylabel('Coherence')
%% Calculate impulse response of optimum filter
[w_opt] = fftshift(ifft(([W_opt; conj(W_opt([1,end:-1:2]))])));     % calcuate impulse response of the optimal control filter
t = -(length(w_opt)/2-1)/Fs:1/Fs:length(w_opt)/2/Fs;                % define time vector including negative time
figure()
plot(t,(w_opt))
xlim([min(t), max(t)])
xlabel('Time, sec')
ylabel('w_{opt}')

%%  <<<<<<<<<<<<<<       PART 4  FxLMS Algorithm    >>>>>>>>>>>>>>>>>>>>>>>

I_w = 128;      % control filter length - you can adjust this to see the effect of control filter length
alpha = 1e-2;   % convergence gain - you can adjust this to see the effect of convergence gain

% Pre-allocate Arrays
w = zeros(I_w, length(d));  % pre-define the control filter - this is a vector of (I_w x 1) at each time step. In order to plot the convergence, I have stored this at each time step
u = zeros(size(d));         % pre-define control signal vector - this is the control signal over time
e = zeros(size(d));         % pre-define error signal vector - this is error signal over time
r = zeros(size(d));         % pre-define filtered reference signal vector - this is the filtered-reference signal over time

% Run LMS
for n=I_w:length(d) % start at I_w
    u(n) = w(:,n).'*x(n:-1:(n-I_w+1)).';            % control signal at time n - filter the reference (current and past I_w samples) with the control filter.
    e(n) = d(n)+G_FIR*u(n:-1:(n-I+1)).';            % error signal at time n - disturbance at time n, G_FIR is the impulse of the plant, u is the vector of current and I past samples of the control signal 
    r(n) = G_FIR*x(n:-1:(n-I+1)).';                 % filtered reference signal - filter reference via plant - the reference vector here contains the current and past I samples
    w(:,n+1)=w(:,n)-alpha*r(n:-1:(n-I_w+1)).'*e(n); % control filter update - filter reference vector here contains the current and past I_w samples
end

% Plot Convergence

figure()
subplot(3,1,1)
plot(t, w(:,1:length(t)))
ylabel('Filter coeffcients')
subplot(3,1,2)
plot(t, u(1:length(t)))
ylabel('Control signal')
subplot(3,1,3)
plot(t, d(1:length(t)),'k',t, e(1:length(t)),'r')
ylabel('Error signal')
xlabel('Time, s')
legend('Disturbance - no control', 'Error - control')

% Plot IR of Converged Filter (Final set)
figure()
plot(0:1/Fs:(I_w-1)/Fs, w(:,end)) 
xlabel('Time, s')
ylabel('w')
axis('tight')

% plot spectrum before and after control
S_dd = pwelch(d(30*Fs:end), win, noverlap, nfft, Fs);
S_ee = pwelch(e(30*Fs:end), win, noverlap, nfft, Fs);
figure() % 
semilogx(f, 10*log10(abs(S_dd)),'k', f, 10*log10(abs(S_ee)), 'r')
xlabel('Frequency, Hz')
xlim([0, max(f)])
ylabel('10log_{10}(S_{ee}), dB')
legend('Disturbance', 'Error')



