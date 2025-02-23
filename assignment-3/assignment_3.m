clear;
close all;

%% A1
N = 100;
bit_seq = (sign(randn(N,4))+1)/2;

%% A2
X = bits_to_PSK_16(bit_seq);
XI = X(:,1)';
XQ = X(:,2)';

%% A3
T = 10^-2;
over = 10;
A = 4;
a = 0.5;
[phi,t] = srrc_pulse(T,over,A,a);

Ts = T/over;
Fs = 1/Ts;
Nf = 2048;
F = -Fs/2:Fs/Nf:Fs/2-Fs/Nf;
XI_delta = Fs*upsample(XI,over);
XQ_delta = Fs*upsample(XQ,over);
t_delta = 0:Ts:N*T-Ts;

XI_t = conv(XI_delta,phi)*Ts;
XQ_t = conv(XQ_delta,phi)*Ts;
t_conv = t_delta(1)+t(1):Ts:t_delta(end)+t(end);

figure;
plot(t_conv,XI_t);
title("Convolution of X_{I,delta}(t) and phi(t)");
xlabel('t(sec)');
ylabel('X_{I}(t)');
grid on;

figure;
plot(t_conv,XQ_t);
title("Convolution of X_{Q,delta}(t) and phi(t)");
xlabel('t(sec)');
ylabel('X_{Q}(t)');
grid on;

XI_f = fftshift(fft(XI_t,Nf));
XI_F = XI_f*Ts;
XQ_f = fftshift(fft(XQ_t,Nf));
XQ_F = XQ_f*Ts;
Ttotal = length(t_conv)*Ts;
PxiF = (abs(XI_F).^2)/Ttotal;
PxqF = (abs(XQ_F).^2)/Ttotal;

figure;
semilogy(F,PxiF);
title("Periodogram of X_{I}(t) using semilogy");
xlabel('F(Hz)');
ylabel('P_{XI}(F)');
grid on;

figure;
semilogy(F,PxqF);
title("Periodogram of X_{Q}(t) using semilogy");
xlabel('F(Hz)');
ylabel('P_{XQ}(F)');
grid on;

%% A4
F0 = 200;
XImod_t = XI_t.*(2*cos(2*pi*F0*t_conv));
XQmod_t = XQ_t.*(-2*sin(2*pi*F0*t_conv));

figure;
plot(t_conv,XImod_t);
title("Modulation X_{I}(t)*(2*cos(2*pi*F0*t))");
xlabel('t(sec)');
ylabel('X_{I,mod}(t)');
grid on;

figure;
plot(t_conv,XQmod_t);
title("Modulation X_{Q}(t)*(-2*sin(2*pi*F0*t))");
xlabel('t(sec)');
ylabel('X_{Q,mod}(t)');
grid on;

XImod_f = fftshift(fft(XImod_t,Nf));
XImod_F = XImod_f*Ts;
XQmod_f = fftshift(fft(XQmod_t,Nf));
XQmod_F = XQmod_f*Ts;
PximodF = (abs(XImod_F).^2)/Ttotal;
PxqmodF = (abs(XQmod_F).^2)/Ttotal;

figure;
semilogy(F,PximodF);
title("Periodogram of X_{I,mod}(t) using semilogy");
xlabel('F(Hz)');
ylabel('P_{XImod}(F)');
grid on;

figure;
semilogy(F,PxqmodF);
title("Periodogram of X_{Q,mod}(t) using semilogy");
xlabel('F(Hz)');
ylabel('P_{XQmod}(F)');
grid on;

%% A5
X_t = XImod_t+XQmod_t;

figure;
plot(t_conv,X_t);
title("Channel Input X(t) = X_{I,mod}(t)+X_{Q,mod}(t)");
xlabel('t(sec)');
ylabel('X(t)');
grid on;

X_f = fftshift(fft(X_t,Nf));
X_F = X_f*Ts;
PxF = (abs(X_F).^2)/Ttotal;

figure;
semilogy(F,PxF);
title("Periodogram of X(t) using semilogy");
xlabel('F(Hz)');
ylabel('P_{X}(F)');
grid on;

%% A6&7
% Add white Gaussian noise W(t) with variance varW to the Channel Output
SNR = 20;
varW = 1/(Ts*(10^(SNR/10)));
W_t = sqrt(varW)*randn(1,length(X_t));
Y_t = X_t+W_t;

%% A8
YIdemod_t = Y_t.*(cos(2*pi*F0*t_conv));
YQdemod_t = Y_t.*(-sin(2*pi*F0*t_conv));

figure;
plot(t_conv,YIdemod_t);
title("Demodulation Υ(t)*(cos(2*pi*F0*t))");
xlabel('t(sec)');
ylabel('Υ_{I,demod}(t)');
grid on;

figure;
plot(t_conv,YQdemod_t);
title("Demodulation Υ(t)*(-sin(2*pi*F0*t))");
xlabel('t(sec)');
ylabel('Υ_{Q,demod}(t)');
grid on;

YIdemod_f = fftshift(fft(YIdemod_t,Nf));
YIdemod_F = YIdemod_f*Ts;
YQdemod_f = fftshift(fft(YQdemod_t,Nf));
YQdemod_F = YQdemod_f*Ts;
PyidemodF = (abs(YIdemod_F).^2)/Ttotal;
PyqdemodF = (abs(YQdemod_F).^2)/Ttotal;

figure;
semilogy(F,PyidemodF);
title("Periodogram of Υ_{I,demod}(t) using semilogy");
xlabel('F(Hz)');
ylabel('P_{YIdemod}(F)');
grid on;

figure;
semilogy(F,PyqdemodF);
title("Periodogram of Υ_{Q,demod}(t) using semilogy");
xlabel('F(Hz)');
ylabel('P_{YQdemod}(F)');
grid on;

%% A9
YI = conv(YIdemod_t,phi)*Ts;
YQ = conv(YQdemod_t,phi)*Ts;
t_conv = t_conv(1)+t(1):Ts:t_conv(end)+t(end);

figure;
plot(t_conv,YI);
title("Convolution of Υ_{I,demod}(t) and phi(t)");
xlabel('t(sec)');
ylabel('Y_{I}(t)');
grid on;

figure;
plot(t_conv,YQ);
title("Convolution of Υ_{Q,demod}(t) and phi(t)");
xlabel('t(sec)');
ylabel('Y_{Q}(t)');
grid on;

YI_f = fftshift(fft(YI,Nf));
YI_F = YI_f*Ts;
YQ_f = fftshift(fft(YQ,Nf));
YQ_F = YQ_f*Ts;
PyiF = (abs(YI_F).^2)/Ttotal;
PyqF = (abs(YQ_F).^2)/Ttotal;

figure;
semilogy(F,PyiF);
title("Periodogram of Y_{I}(t) using semilogy");
xlabel('F(Hz)');
ylabel('P_{YI}(F)');
grid on;

figure;
semilogy(F,PyqF);
title("Periodogram of Y_{Q}(t) using semilogy");
xlabel('F(Hz)');
ylabel('P_{YQ}(F)');
grid on;

%% A10
YI_sampling = YI((2*A*T/Ts)+1:over:length(YI)-(2*A*T/Ts));
YQ_sampling = YQ((2*A*T/Ts)+1:over:length(YQ)-(2*A*T/Ts));

Y_sampling = zeros(N,2);
for i = 1:N
    Y_sampling(i,1) = YI_sampling(i);
    Y_sampling(i,2) = YQ_sampling(i);
end

scatterplot(Y_sampling);
title("Sampling of Y_{I}(t) and Y_{Q}(t)");
grid on;

%% A11
[est_X,est_bit_seq] = detect_PSK_16(Y_sampling);

%% A12
num_of_symbol_errors = symbol_errors(est_X,X);

%% A13
num_of_bit_errors = bit_errors(est_bit_seq,bit_seq);
