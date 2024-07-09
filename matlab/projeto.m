%% PROJETO CONTROLADOR RESSONANTE PARA VSC1 (FORMADOR DE REDE) %%
clc
clear all
close all
s = tf('s');

%% TEMPOS E FREQUÊNCIAS
% do sistema
f_rede = 60;
t_rede = 1/f_rede;
f = 20e3;
ts = 1/f;
omega0 = 2*pi*f_rede;
% vetor de tempo
t_final = 10*t_rede;
t = 0:ts:t_final;
% sinal de referência
sin_w60 = sin(omega0*t);

%% PARÂMETROS DO SISTEMA
% indutancia do filtro
Rf = 0.1;
% resistencia parasita
Lf = 1e-3;
% capacitancia do filtro
Cf = 44e-6;
% resistencia de amortecimento passivo
Rd = 1;
% impedancia da linha
Rl = 10e-3;
Ll = 10e-3;
% impedancia da carga + linha (S = 12000 + 9000i)
L0 = 0.0051 + Ll;
R0 = 2.5813 + Rl;
% L0 = Ll;
% R0 = Rl;

%% MODELAGEM DO SISTEMA
% X = [I_{Lf} V_{c} I_{L0}]'
A = [(-Rf-Rd)/Lf -1/Lf Rd/Lf; 1/Cf 0 -1/Cf; Rd/L0 1/L0 (-R0-Rd)/L0];
B = [1/Lf; 0; 0];
C_Ilf = [1 0 0];
C_Vc = [0 1 0];
C_I0 = [0 0 1];
D = 0;

%% EXTRACAO DE FUNCOES DE TRANSFERENCIA: Continuo
sys_ilf = ss(A,B,C_Ilf,D);
[num, den] = ss2tf(A,B,C_Ilf,D);
gp_ilf = tf(num, den);

sys_vc = ss(A,B,C_Vc,D);
[num, den] = ss2tf(A,B,C_Vc,D);
gp_vc = tf(num, den);

sys_i0 = ss(A,B,C_I0,D);
[num, den] = ss2tf(A,B,C_I0,D);
gp_i0 = tf(num, den);

%% DISCRETIZACAO DO SISTEMA
sysd = c2d(sys_ilf, ts, 'zoh');
[Ad, Bd, Cd, Dd] = ssdata(sysd);

A = [Ad, Bd; 0 0 0 0];
B = [0; 0; 0; 1];
C_Ilf = [1 0 0 0];
C_Vc =  [0 1 0 0];
C_I0 =  [0 0 1 0];
D = 0;
sysd = ss(A, B, C_Ilf, D, ts);

%% PLOTAGEM: Contínuo vs Discreto (+ Atraso)
figure(1)
step(sys_ilf, sysd);
grid on;
legend('continuo', 'discreto');

%% EXTRACAO DE FUNCOES DE TRANSFERENCIA: Discreto
[num, den] = ss2tf(A,B,C_Ilf,D);
gp_ilfd = tf(num, den, ts);

[num, den] = ss2tf(A,B,C_Vc,D);
gp_vcd = tf(num, den, ts);

[num, den] = ss2tf(A,B,C_I0,D);
gp_i0d = tf(num, den, ts);

%% PROJETO: Controlador de Corrente
Q = [1000 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1];
R = 1;
K = dlqr(A, B, Q, R)

%% Sistema Realimentado
A = A-B*K;
B = [0; 0; 0; 1*14.48];
sysd_mf = ss(A, B, C_Ilf, D, ts);

figure(2)
bode(sysd, sysd_mf)
grid on;
legend('nominal', 'realimentado');

%% PROJETO: Controlador Ressonante para Tensï¿½o

a = 0.1;
R = [0 1; -exp(-2*a*ts) 2*exp(-a*ts)*cos(omega0*ts)];
T = [0; 1];

%% PROJETO: Controlador de Tensao
Q = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1];
R = 10000;
K = dlqr(A, B, Q, R)