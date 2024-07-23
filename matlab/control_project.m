%% RESONANT CONTROLLER DESIGN FOR VSC1 (GRID FORMING)
clc;
clear all;
close all;
s = tf('s');

%% TIMES AND FREQUENCIES
% system
f_grid = 60;
t_grid = 1/f_grid;
f = 20e3;
ts = 1/f;
omega0 = 2*pi*f_grid;
% time vector
t_final = 10*t_grid;
t = 0:ts:t_final;
% reference signal
sin_w60 = sin(omega0*t);
% z
z = tf('z', ts);

%% SYSTEM PARAMETERS
% filter inductance
Rf = 0.1;
% parasitic resistance
Lf = 1e-3;
% filter capacitance
Cf = 44e-6;
% passive damping resistance
Rd = 0.5;
% load impedance (S = 12000 + 9000i)
L0 = 0.0051;
R0 = 2.5813;

%% SYSTEM MODELING
% X = [I_{Lf} V_{c} I_{L0}]'
A = [(-Rf-Rd)/Lf -1/Lf Rd/Lf; 1/Cf 0 -1/Cf; Rd/L0 1/L0 (-R0-Rd)/L0];
B = [1/Lf; 0; 0];
C = [0 1 0];
D = 0;
sys = ss(A, B, C, D);

%% EXTRACTION OF TRANSFER FUNCTIONS: Continuous
gc = C * inv((s*eye(3,3) - A)) * B;

%% SYSTEM DISCRETIZATION
Ad = expm(A*ts);
Bd = (Ad - eye(3,3)) * inv(A) * B;
A = [Ad, Bd; 0 0 0 0];
B = [0; 0; 0; 1];
C = [0 1 0 0];
D = 0;
sysd = ss(A, B, C, D, ts);

%% EXTRACTION OF TRANSFER FUNCTIONS: Discrete
gd = C * inv((z*eye(4,4) - A)) * B;
% extraction of simplified function
[N, D] = tfdata(gd);
beta2 = N{1}(3);
beta1 = N{1}(4);
beta0 = N{1}(5);
alpha3 = D{1}(2);
alpha2 = D{1}(3);
alpha1 = D{1}(4);
num = (beta1 + 2*beta2)*z + beta0 - beta2;
den = z^4 + alpha3 * z^3 + alpha2 * z^2 + alpha1*z;
g0 = num/den;

%% PLANT DIVISION INTO TWO DYNAMICS: g0 and muDelta
z1 = 9.750106713514798e-01;
z2 = -9.884000680209557e-01;
k2 = beta2;
k1 = beta1 - 2*k2*z2;
k0 = beta0 - k1*z1 - k2 * (z2^2);
g0 = (k1 * (z+z1) + k0)/den;
mu = (k2 * (z+z2)^2)/den;

%% REFERENCE MODEL
p = 0.3;
num = (1-p)^3;
den = (z-p)^3;
wmd = tf(num/den, 1/20000);

%% PLOTTING: bode of g, g0, muDelta
figure(1);
bode = bodeplot(gd, g0, mu, wmd);
title('Bode Diagram');
setoptions(bode, 'FreqScale', 'log', 'FreqUnits', 'Hz', 'Grid', 'on');
sub = legend('$G(z)$', '$G\_0(z)$', '$\mu\Delta\_a(z)$', '$W\_m(z)$');
set(sub, 'Interpreter', 'latex');
set(sub, 'Location', 'southwest');
xlim([1 1e4]);
xlabel('Frequency');
ylabel('Phase');

%% DESIGN: virtual impedance
% electrical parameters
lv = 2e-3;
rv = 1e-4;
% parameters of zeros
wp = 2*pi*1500;
zeta = 1;
zo = s*lv;
% IDEAL:
zv0 = rv + s*lv;
% REALIZABLE
zv = ((rv + s*lv)*wp^2) / (s^2 + 2*zeta*wp*s + wp^2);

%% DESIGN: low pass filter
wc = 2*pi*800;
lpf = wc/(s+wc);

%% PLOTTING: virtual impedance
figure(2);
bode = bodeplot(zv0, zv, lpf*zv);
setoptions(bode, 'FreqScale', 'log', 'FreqUnits', 'Hz', 'Grid', 'on');
title('Bode Diagram');
sub = legend('$Z\_{v\_0}$', '$Z\_{v}$', '$LPF\_{800Hz} \cdot Z\_{v}$');
set(sub, 'Interpreter', 'latex');
set(sub, 'Location', 'southwest');
grid on;
xlim([10 10e3]);