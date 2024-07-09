%% PROJETO CONTROLADOR RESSONANTE PARA VSC1 (FORMADOR DE REDE) %%
clc
clear all
close all
s = tf('s');
z = tf('z', 1/20000);
format long e;

%% PROJETO: Modelo de Referência
settling_time = 2e-3;
tau = settling_time/4;
pole_1 = 1/tau
pole_2 = 5 * pole_1
pole_3 = 5 * pole_1

num = pole_1 * pole_2 * pole_3;
den = (s + pole_1) * (s + pole_2) * (s + pole_3);
wm = num/den;

%% STEP: Modelo de Referência
figure(1)
step(wm);

%% DISCRETIZAÇÃO: Modelo de Referência
p = 0.25;
num = (1-p)^3;
den = (z-p)^3;
wmd = tf(num/den, 1/20000);

%wmd = c2d(wm, 1/20000, 'tustin');
[num, den] = tfdata(wmd, 'v');

% Extract and print numerator (alpha) coefficients
for i = 1:length(num)
    fprintf('ALPHA%d = %+1.15e;\n', i-1, num(i));
end

% Extract and print denominator (beta) coefficients
for i = 2:length(den) % Skip the first coefficient because it's always 1 for tf
    fprintf('BETA%d = %+1.15e;\n', i-1, den(i));
end

figure(2);
bode(wmd);
xlim([1 10e3]);

figure(3);
step(wmd);

pho = 1/db2mag(-3.7e-3)
