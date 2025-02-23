clear;
close all;

%% A1
T = 10^-2;
over = 10;
A = 4;
a = 0.5;
[phi,t] = srrc_pulse(T,over,A,a);
figure;
plot(t,phi);
title("Square root raised cosine pulse");
xlim([-4*T 4*T]);
xlabel('t(sec)');
ylabel('phi(t)');
grid on;

Ts = T/over;
Fs = 1/Ts;
Nf = 2048;
F = -Fs/2:Fs/Nf:Fs/2-Fs/Nf;
PHI_f = fftshift(fft(phi,Nf));
PHI_F = PHI_f*Ts;
figure;
semilogy(F,abs(PHI_F).^2);
title("Energy Spectral Density of SRRC Pulse using semilogy");
xlabel('F(Hz)');
ylabel('|PHI(F)|^2');
grid on;

%% A2
N = 100;
b = (sign(randn(N,1))+1)/2; % Random sequence of N bits
n = 0:N-1;
X = bits_to_2PAM(b);
figure;
hold on;
stem(n,b);
stem(n,X);
hold off;
legend('b (bits)','X (2-PAM symbols)');

X_delta = 1/Ts*upsample(X,over);
t_delta = 0:Ts:N*T-Ts;
figure;
stem(t_delta,X_delta);
xlabel('t(sec)');
ylabel('X_{delta}(t)');

X_t = conv(X_delta,phi)*Ts;
t_conv = t_delta(1)+t(1):Ts:t_delta(end)+t(end);
figure;
plot(t_conv,X_t);
title("Convolution of X_{delta}(t) and phi(t)");
xlabel('t(sec)');
ylabel('X(t)');
grid on;

%% A3
X_f = fftshift(fft(X_t,Nf));
X_F = X_f*Ts;
Ttotal = length(t_conv)*Ts;
PxF = (abs(X_F).^2)/Ttotal;
figure;
plot(F,PxF);
title("Periodogram of X(t) using plot");
xlabel('F(Hz)');
ylabel('P_{X}(F)');
grid on;
figure;
semilogy(F,PxF);
title("Periodogram of X(t) using semilogy");
xlabel('F(Hz)');
ylabel('P_{X}(F)');
grid on;

% N = 50;
% K = 250;
N = 100;
K = 500;
PxF_K = zeros(K,Nf);
for i = 1:K
    b = (sign(randn(N,1))+1)/2;
    X = bits_to_2PAM(b);
    X_delta = 1/Ts*upsample(X,over);
    X_t = conv(X_delta,phi)*Ts;
    X_F = fftshift(fft(X_t,Nf))*Ts;
    PxF_K(i,:) = (abs(X_F).^2)/Ttotal;
end
SxF_exp = mean(PxF_K,1);
SxF_theor = (var(X)/T).*(abs(PHI_F).^2);
figure;
semilogy(F,SxF_exp);
hold on;
semilogy(F,SxF_theor);
hold off;
legend('Experimental','Theoretical');
title("Power Spectral Density of X(t) using semilogy");
xlabel('F(Hz)');
ylabel('S_{X}(F)');
grid on;

%% A4
b = (sign(randn(N/2,2))+1)/2;
n = 0:N/2-1;
X = bits_to_4PAM(b);
figure;
stem(n,X);
title("X (4-PAM symbols)");

X_delta = 1/Ts*upsample(X,over);
t_delta = 0:Ts:(N/2)*T-Ts;
X_t = conv(X_delta,phi)*Ts;
t_conv = t_delta(1)+t(1):Ts:t_delta(end)+t(end);
figure;
plot(t_conv,X_t);
title("Convolution of X_{delta}(t) and phi(t)");
xlabel('t(sec)');
ylabel('X(t)');
grid on;

X_f = fftshift(fft(X_t,Nf));
X_F = X_f*Ts;
Ttotal = length(t_conv)*Ts;
PxF = (abs(X_F).^2)/Ttotal;
figure;
semilogy(F,PxF);
title("Periodogram of X(t) using semilogy");
xlabel('F(Hz)');
ylabel('P_{X}(F)');
grid on;

N = 100;
K = 500;
PxF_K = zeros(K,Nf);
for i = 1:K
    b = (sign(randn(N/2,2))+1)/2;
    X = bits_to_4PAM(b);
    X_delta = 1/Ts*upsample(X,over);
    X_t = conv(X_delta,phi)*Ts;
    X_F = fftshift(fft(X_t,Nf))*Ts;
    PxF_K(i,:) = (abs(X_F).^2)/Ttotal;
end
SxF_exp = mean(PxF_K,1);
SxF_theor = (var(X)/T).*(abs(PHI_F).^2);
figure;
semilogy(F,SxF_exp);
hold on;
semilogy(F,SxF_theor);
hold off;
legend('Experimental','Theoretical');
title("Power Spectral Density of X(t) using semilogy");
xlabel('F(Hz)');
ylabel('S_{X}(F)');
grid on;

%% A5
T = 2*10^-2;
over = 20;
A = 4;
a = 0.5;
[phi,t] = srrc_pulse(T,over,A,a);
Ts = T/over;
Fs = 1/Ts;
Nf = 2048;
F = -Fs/2:Fs/Nf:Fs/2-Fs/Nf;
PHI_f = fftshift(fft(phi,Nf));
PHI_F = PHI_f*Ts;
N = 100;
b = (sign(randn(N,1))+1)/2;
X = bits_to_2PAM(b);
X_delta = 1/Ts*upsample(X,over);
t_delta = 0:Ts:N*T-Ts;
X_t = conv(X_delta,phi)*Ts;
t_conv = t_delta(1)+t(1):Ts:t_delta(end)+t(end);
X_f = fftshift(fft(X_t,Nf));
X_F = X_f*Ts;
Ttotal = length(t_conv)*Ts;
PxF = (abs(X_F).^2)/Ttotal;
figure;
plot(F,PxF);
title("Periodogram of X(t) using plot");
xlabel('F(Hz)');
ylabel('P_{X}(F)');
grid on;
figure;
semilogy(F,PxF);
title("Periodogram of X(t) using semilogy");
xlabel('F(Hz)');
ylabel('P_{X}(F)');
grid on;
% N = 50;
% K = 250;
N = 100;
K = 500;
PxF_K = zeros(K,Nf);
for i = 1:K
    b = (sign(randn(N,1))+1)/2;
    X = bits_to_2PAM(b);
    X_delta = 1/Ts*upsample(X,over);
    X_t = conv(X_delta,phi)*Ts;
    X_F = fftshift(fft(X_t,Nf))*Ts;
    PxF_K(i,:) = (abs(X_F).^2)/Ttotal;
end
SxF_exp = mean(PxF_K,1);
SxF_theor = (var(X)/T).*(abs(PHI_F).^2);
figure;
semilogy(F,SxF_exp);
hold on;
semilogy(F,SxF_theor);
hold off;
legend('Experimental','Theoretical');
title("Power Spectral Density of X(t) using semilogy");
xlabel('F(Hz)');
ylabel('S_{X}(F)');
grid on;

%% B1
dt = 0.00001;
T = 0.001;
Fo = 1000;
t = -5*T:dt:5*T;
k = 5;
N = 0-dt:dt:1-dt;
U = 0:dt:2*pi-dt;
figure; 
hold on;
colors = {'r', 'g', 'b', 'm', 'c'};
labels = cell(1, k);
for i = 1:k
    % Select random values for X and phi
    permuted_indices_N = randperm(length(N));
    random_indices_N = permuted_indices_N(1);
    X = N(random_indices_N);
    permuted_indices_U = randperm(length(U));
    random_indices_U = permuted_indices_U(1);
    phi = U(random_indices_U);
    Y = X*cos(2*pi*Fo*t + phi);
    plot(t, Y, 'color', colors{i});
    labels{i} = sprintf('For X=%.3f, {{\\phi}} = %.3f', X, phi);
end
hold off;
xlabel('Time (sec)');
ylabel('Amplitude ( X(t) )');
title('Stochastic processes" Y(t) = X(t)*cos(2*{\pi}*F0*t + {\phi})" ');
legend(labels);
