%% RESONANT CONTROLLER DESIGN FOR VSC1 (GRID FORMING) %%
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
C = [1 0 0];
D = 0;
sys = ss(A, B, C, D);

%% LOADING PSIM .MAT FILE
load('model_validation.mat');

%% STEP RESPONSE
tsc = 1e-7;
t = 0:tsc:20e-3;
[y, t] = step(sys, t);
plot(t, y, '.');
hold on;
plot(Time, I1);
legend('MATLAB', 'PSIM');