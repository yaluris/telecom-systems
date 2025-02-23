clear;
close all;

%% TH1
T = 1;
t = -2:0.01:2;
R = zeros(1,length(t));
for i = 1:length(t)
    if (abs(t(i)) <= T)
        R(i) = 1-abs(t(i))/T;
    else
        R(i) = 0;
    end
end
figure;
plot(t,R);
title("R_{φφ}(τ)");
grid on;

%% TH3
T = 1;
t = -2:0.01:2;
R = zeros(1,length(t));
for i = 1:length(t)
    if (t(i) > T/2) && (t(i) <= T)
        R(i) = -1+t(i)/T;
    elseif (t(i) >= -T) && (t(i) < -T/2)
        R(i) = -1-t(i)/T;
    elseif (t(i) > 0) && (t(i) <= T/2)
        R(i) = 1-3*t(i)/T;
    elseif (t(i) >= -T/2) && (t(i) <= 0)
        R(i) = 1+3*t(i)/T;
    else
        R(i) = 0;
    end
end
figure;
plot(t,R);
title("R_{φφ}(τ)");
grid on;

%% A1
T = 10^-3;
over = 10;
A = 4;
a1 = 0;
a2 = 0.5;
a3 = 1;
[phi1,t] = srrc_pulse(T,over,A,a1);
[phi2,t] = srrc_pulse(T,over,A,a2);
[phi3,t] = srrc_pulse(T,over,A,a3);
figure;
plot(t,phi1);
hold on;
plot(t,phi2);
plot(t,phi3);
hold off;
legend('a = 0','a = 0.5','a = 1');
xlim([-4*T 4*T]);
title("Square Root Raised Cosine Pulses");
xlabel('t(sec)');
ylabel('phi(t)');
grid on;

%% A2
Ts = T/over;
Fs = 1/Ts;
N = 2048;
F = [-Fs/2:Fs/N:Fs/2-Fs/N];
PHI1_F = fftshift(fft(phi1,N)*Ts);
PHI2_F = fftshift(fft(phi2,N)*Ts);
PHI3_F = fftshift(fft(phi3,N)*Ts);
figure;
plot(F,abs(PHI1_F).^2);
hold on;
plot(F,abs(PHI2_F).^2);
plot(F,abs(PHI3_F).^2);
hold off;
legend('a = 0','a = 0.5','a = 1');
title("Energy Spectral Density of SRRC Pulses using plot");
xlabel('F(Hz)');
ylabel('|PHI(F)|^2');
grid on;

figure;
semilogy(F,abs(PHI1_F).^2);
hold on;
semilogy(F,abs(PHI2_F).^2);
semilogy(F,abs(PHI3_F).^2);
hold off;
legend('a = 0','a = 0.5','a = 1');
title("Energy Spectral Density of SRRC Pulses using semilogy");
xlabel('F(Hz)');
ylabel('|PHI(F)|^2');
grid on;

%% A3
% c = (T/10^3)*ones(1,N);
c = (T/10^5)*ones(1,N);
figure;
semilogy(F,abs(PHI1_F).^2);
hold on;
semilogy(F,abs(PHI2_F).^2);
semilogy(F,abs(PHI3_F).^2);
semilogy(F,c);
hold off;
legend('a = 0','a = 0.5','a = 1');
title("Energy Spectral Density of SRRC Pulses using semilogy");
xlabel('F(Hz)');
ylabel('|PHI(F)|^2');
grid on;

%% C1
T = 10^-3;
over = 10;
Ts = T/over;
A = 4;
a = 0.5;
N = 50;
b = (sign(randn(N,1))+1)/2;
n = 0:N-1;
figure;
stem(n,b);
title("Random sequence of N bits");

%% C2
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

[phi,t] = srrc_pulse(T,over,A,a);
figure;
plot(t,phi);
title("SRRC Pulse");
xlim([-4*T 4*T]);
xlabel('t(sec)');
ylabel('phi(t)');
grid on;

X_t = conv(X_delta,phi)*Ts;
t_conv = t_delta(1)+t(1):Ts:t_delta(end)+t(end);
figure;
plot(t_conv,X_t);
title("Convolution of X_{delta}(t) and phi(t)");
xlabel('t(sec)');
ylabel('X(t)');
grid on;

Z = conv(X_t,phi)*Ts;
t_z = t_conv(1)+t(1):Ts:t_conv(end)+t(end);
figure;
plot(t_z,Z);
hold on;
stem([0:N-1]*T,X);
hold off;
title("Convolution of X(t) and phi(-t)");
legend('Z(t)','Z(kT), k = 0,...,N−1');
xlabel('t(sec)');
grid on;

%% B1
T = 10^-3;
over = 10;
A = 4;
a1 = 0;
a2 = 0.5;
a3 = 1;

% Calculating SRRC pulses
[phi1, t] = srrc_pulse(T, over, A, a1);
[phi2, t] = srrc_pulse(T, over, A, a2);
[phi3, t] = srrc_pulse(T, over, A, a3);

k0 = 0;
k1 = 1;
k2 = 2;
k3 = 3;

% Time shift
tShift0 = t + (k0 * T);   % phi(t)
tShift1 = t + (k1 * T);   % phi(t-T)
tShift2 = t + (k2 * T);   % phi(t-2T)
tShift3 = t + (k3 * T);   % phi(t-3T)

% Compute p1T, p2T, p3T
p1T = 4*a3/(pi*sqrt(T)) * phi1 .* cos((1+a3)*pi*(t-T)/T) ./ (1-(4*a3*(t-T)./T).^2);
p2T = 4*a3/(pi*sqrt(T)) * phi1 .* cos((1+a3)*pi*(t-2*T)/T) ./ (1-(4*a3*(t-2*T)./T).^2);
num = cos((1+a3)*pi*(t-k3*T)/T) + sin((1-a3)*pi*(t-k3*T)/T) ./ (4*a3*(t-k3*T)/T);
denom = 1-(4*a3*(t-k3*T)./T).^2;
p3T = 4*a3/(pi*sqrt(T)) * num ./ denom;

% Plot for k=0
figure
subplot(2,1,1)
hold on
plot(t, phi3, 'linewidth', 1.5)   %
plot(t, phi3, 'linewidth', 1.5)   %
legend('\phi(t-0*T)');  %
set(gca, 'fontsize', 15)
title('For a = 1, k = 0');      %
grid on
ylabel('Shifted SRRC pulses ');
xlabel('Time in sec');
pro = phi3.*phi3;   %
subplot(2,1,2)
plot(t, pro, 'linewidth', 1.5)
set(gca, 'fontsize', 15)
title('For a = 1, k = 0');
grid on
ylabel('Product of shifted  SRRC pulses ');
xlabel('Time in sec');

% Numerical integration of data using the trapezoidal method. 
int0 =  trapz(t,pro)

%%%%%%For k = 1%%%%%%%%%%%%%%%%%

figure
subplot(2,1,1)
hold on
plot(t, phi3, 'linewidth', 1.5)   %
plot(t, p1T, 'linewidth', 1.5)    %
legend('\phi(t-1*T)');
set(gca, 'fontsize', 15)
title('For a = 1, k = 1');
grid on
ylabel('Shifted SRRC pulses ');
xlabel('Time in sec');

pro = phi3.*p1T;    %
subplot(2,1,2)
plot(t, pro, 'linewidth', 1.5)
set(gca, 'fontsize', 15)
title('For a = 1, k = 1');
grid on
ylabel('Product of shifted  SRRC pulses ');
xlabel('Time in sec');

int1 =  trapz(t,pro)
%%%%%%For k = 2%%%%%%%%%%%%%%%%%

figure
subplot(2,1,1)
hold on
plot(t, phi3, 'linewidth', 1.5)   %
plot(t, p2T, 'linewidth', 1.5)
legend('\phi(t-2*T)');
set(gca, 'fontsize', 15)
title('For a = 1, k = 2');
grid on
ylabel('Shifted SRRC pulses ');
xlabel('Time in sec');

pro = phi3.*p2T;    %
subplot(2,1,2)
plot(t, pro, 'linewidth', 1.5)
set(gca, 'fontsize', 15)
title('For a = 1, k = 2');
grid on
ylabel('Product of shifted  SRRC pulses ');
xlabel('Time in sec');

int2 =  trapz(t,pro)
%%%%%%For k = 3%%%%%%%%%%%%%%%%%

figure
subplot(2,1,1)
hold on
plot(t, phi3, 'linewidth', 1.5)   %
plot(t, p3T, 'linewidth', 1.5)
legend('\phi(t-3*T)');
set(gca, 'fontsize', 15)
title('For a = 1, k = 3');
grid on
ylabel('Shifted SRRC pulses ');
xlabel('Time in sec');

pro = phi3.*p3T;    %
subplot(2,1,2)
plot(t, pro, 'linewidth', 1.5)
set(gca, 'fontsize', 15)
title('For a = 1, k = 3');
grid on
ylabel('Product of shifted  SRRC pulses ');
xlabel('Time in sec');

int3 =  trapz(t,pro)
