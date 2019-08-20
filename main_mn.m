
clear all; clc;

Tmax = 0.2;
dT = 10^(-4);
ovhd = 5; %mV

nInt = 2;
%tvec = 1:maxIter;

sgtitle('Mihalas-Niebur Neuron Features')

%% A - TONIC SPIKING %% 
Tmax = 0.2;
tvec = 0:dT:Tmax;
a = 0; % s^-1
A = [0, 0]; % V/s
% Ie = 1.5; % V/s
Ie = 1.5*ones(1, length(tvec));
[Iint, Vm, Theta] = MN_Calc(a, A, Ie, Tmax, nInt);

% Iints(1, :, :) = Iint;
% Vms(1, :) = Vm;
% Thetas(1, :) = Theta;

subplot(4, 5, 1)
yyaxis left; plot(tvec, Vm) 
ylabel('memrane potential (mV)')
ylim([1.1*min(Vm), max(Vm)+ovhd])
yyaxis right; plot(tvec, Ie) 
title('A: Tonic Spiking')


%% B - CLASS 1 %%
Tmax = 0.5;
tvec = 0:dT:Tmax;
a = 0; % s^-1
A = [0, 0]; % V/s
% Ie = 1.5; % V/s
Ie = (1+10^(-6))*ones(1, length(tvec));
[Iint, Vm, Theta] = MN_Calc(a, A, Ie, Tmax, nInt);

% Iints(2, :, :) = Iint;
% Vms(2, :) = Vm;
% Thetas(2, :) = Theta;

subplot(4, 5, 2)
yyaxis left; plot(tvec, Vm) 
ylim([1.1*min(Vm), max(Vm)+ovhd])
yyaxis right; plot(tvec, Ie) 
title('B: Class 1')

%% C - SPIKE FREQUENCY ADAPTATION %%
Tmax = 0.2;
tvec = 0:dT:Tmax;
a = 5; % s^-1
A = [0, 0]; % V/s
% Ie = 1.5; % V/s
Ie = 2*ones(1, length(tvec));
[Iint, Vm, Theta] = MN_Calc(a, A, Ie, Tmax, nInt);

% Iints(3, :) = Iint;
% Vms(3, :) = Vm;
% Thetas(3, :) = Theta;

subplot(4, 5, 3)
yyaxis left; plot(tvec, Vm) 
ylim([1.1*min(Vm), max(Vm)+ovhd])
yyaxis right; plot(tvec, Ie) 
title('C: Spike Frequency Adaptation')

%% D - PHASIC SPIKING %%
Tmax = 0.5;
tvec = 0:dT:Tmax;
a = 5; % s^-1
A = [0, 0]; % V/s
% Ie = 1.5; % V/s
Ie = 1.5*ones(1, length(tvec));
[Iint, Vm, Theta] = MN_Calc(a, A, Ie, Tmax, nInt);

% Iints(4, :) = Iint;
% Vms(4, :) = Vm;
% Thetas(4, :) = Theta;

subplot(4, 5, 4)
yyaxis left; plot(tvec, Vm) 
ylim([1.1*min(Vm), max(Vm)+ovhd])
yyaxis right; plot(tvec, Ie) 
title('D: Phasic Spiking')

%% E - ACCOMODATION %% 
Tmax = 1;
tvec = 0:dT:Tmax;
LT = length(tvec);
a = 5; % s^-1
A = [0, 0]; % V/s
% Ie = 1.5; % V/s
Ie = ones(1, length(tvec));

Ie(1:LT/6) = 1.5;
Ie(LT/6+1:2*LT/6) = 0;
Ie(2*LT/6+1:3*LT/6) = 0.5;
Ie(3*LT/6+1:4*LT/6) = 1;
Ie(4*LT/6+1:5*LT/6) = 1.5;
Ie(5*LT/6+1:end) = 0;

[Iint, Vm, Theta] = MN_Calc(a, A, Ie, Tmax, nInt);

% Iints(5, :) = Iint;
% Vms(5, :) = Vm;
% Thetas(5, :) = Theta;

subplot(4, 5, 5)
yyaxis left; plot(tvec, Vm)
ylim([1.1*min(Vm), max(Vm)+ovhd])
yyaxis right; plot(tvec, Ie) 
ylim([-1, 2.5])
title('E: Accomodation')
ylabel('external current (V/s)')

%% F - THRESHOLD VARIABILITY %% 
Tmax = 0.4;
tvec = 0:dT:Tmax;
LT = length(tvec);
a = 5; % s^-1
A = [0, 0]; % V/s
% Ie = 1.5; % V/s
Ie = ones(1, length(tvec));

Ie(1:LT/6) = 1.5;
Ie(LT/6+1:2*LT/6) = 0;
Ie(2*LT/6+1:3*LT/6) = -1.5;
Ie(3*LT/6+1:4*LT/6) = 0;
Ie(4*LT/6+1:5*LT/6) = 1.5;
Ie(5*LT/6+1:end) = 0;

[Iint, Vm, Theta] = MN_Calc(a, A, Ie, Tmax, nInt);

% Iints(6, :) = Iint;
% Vms(6, :) = Vm;
% Thetas(6, :) = Theta;

subplot(4, 5, 6)
yyaxis left; plot(tvec, Vm) 
ylim([1.1*min(Vm), max(Vm)+ovhd])
ylabel('membrane potential (mV)')
yyaxis right; plot(tvec, Ie) 
ylim([-2.5, 2.5])
title('F: Threshold Variability')

%% G - REBOUND SPIKE %% 
Tmax = 1;
tvec = 0:dT:Tmax;
LT = length(tvec);
a = 5; % s^-1
A = [0, 0]; % V/s
% Ie = 1.5; % V/s
Ie = ones(1, length(tvec));

Ie(1:LT/3) = 0;
Ie(LT/3+1:2*LT/3) = -3.5;
Ie(2*LT/3+1:end) = 0;

[Iint, Vm, Theta] = MN_Calc(a, A, Ie, Tmax, nInt);

% Iints(7, :) = Iint;
% Vms(7, :) = Vm;
% Thetas(7, :) = Theta;

subplot(4, 5, 7)
yyaxis left; plot(tvec, Vm) 
ylim([1.1*min(Vm), max(Vm)+ovhd])
yyaxis right; plot(tvec, Ie) 
ylim([-4.5, 1])
title('G: Rebound Spike')


%% H - CLASS 2 %% 
Tmax = 0.3;
tvec = 0:dT:Tmax;
LT = length(tvec);
a = 5; % s^-1
A = [0, 0]; % V/s
% Ie = 1.5; % V/s
Ie = 2*(1+10^(-6))*ones(1, length(tvec));
[Iint, Vm, Theta] = MN_Calc(a, A, Ie, Tmax, nInt);

subplot(4, 5, 8)
yyaxis left; plot(tvec, Vm) 
ylim([1.1*min(Vm), max(Vm)+ovhd])
yyaxis right; plot(tvec, Ie) 
title('H: Class 2')


%% I - INTEGRATOR %% 
Tmax = 0.4;
tvec = 0:dT:Tmax;
LT = length(tvec);
a = 5; % s^-1
A = [0, 0]; % V/s
% Ie = 1.5; % V/s
Ie = ones(1, length(tvec));
Ie(1:LT/8) = 1.5;
Ie(LT/8+1:2*LT/8) = 0;
Ie(2*LT/8+1:3*LT/8) = 1.5;
Ie(3*LT/8+1:4*LT/8) = 0;
Ie(4*LT/8+1:5*LT/8) = 1.5;
Ie(5*LT/8+1:6*LT/8) = 0;
Ie(6*LT/8+1:7*LT/8) = 0;
Ie(7*LT/8+1:end) = 0;

[Iint, Vm, Theta] = MN_Calc(a, A, Ie, Tmax, nInt);

subplot(4, 5, 9)
yyaxis left; plot(tvec, Vm) 
ylim([1.1*min(Vm), max(Vm)+ovhd])
yyaxis right; plot(tvec, Ie) 
ylim([-1, 2.5])
ylabel([-1, 2.5])
title('I: Integrator')

%% J - INPUT BISTABILITY %% 
Tmax = 1;
tvec = 0:dT:Tmax;
LT = length(tvec);
a = 5; % s^-1
A = [0, 0]; % V/s
% Ie = 1.5; % V/s
Ie = ones(1, length(tvec));
Ie(1:LT/4) = 1.5;
Ie(LT/4+1:2*LT/4) = 1.7;
Ie(2*LT/4+1:3*LT/4) = 1.5;
Ie(3*LT/4+1:end) = 1.7;

[Iint, Vm, Theta] = MN_Calc(a, A, Ie, Tmax, nInt);

subplot(4, 5, 10)
yyaxis left; plot(tvec, Vm) 
ylim([1.1*min(Vm), max(Vm)+ovhd])
yyaxis right; plot(tvec, Ie) 
ylim([1, 2])
title('J: Input Bistability')
ylabel('external current (V/s)')

%% K - HYPERPOLARIZING SPIKING %%
Tmax = 0.4;
tvec = 0:dT:Tmax;
LT = length(tvec);
a = 30; % s^-1
A = [0, 0]; % V/s
% Ie = 1.5; % V/s
Ie = -1*ones(1, length(tvec));

[Iint, Vm, Theta] = MN_Calc(a, A, Ie, Tmax, nInt);

subplot(4, 5, 11)
yyaxis left; plot(tvec, Vm) 
ylim([1.1*min(Vm), max(Vm)+ovhd])
ylabel('membrane potential (mV)')
yyaxis right; plot(tvec, Ie) 
title('K: Hyperpolarizing Spiking')


%% L - HYPERPOLARIZING BURSTING %%
Tmax = 0.4;
tvec = 0:dT:Tmax;
LT = length(tvec);
a = 30; % s^-1
A = [10, -0.6]; % V/s
% Ie = 1.5; % V/s
Ie = -1*ones(1, length(tvec));

[Iint, Vm, Theta] = MN_Calc(a, A, Ie, Tmax, nInt);

subplot(4, 5, 12)
yyaxis left; plot(tvec, Vm) 
ylim([1.1*min(Vm), max(Vm)+ovhd])
yyaxis right; plot(tvec, Ie) 
title('L: Hyperpolarizing Bursting')

%% M - TONIC BURSTING %%
Tmax = 0.5;
tvec = 0:dT:Tmax;
LT = length(tvec);
a = 5; % s^-1
A = [10, -0.6]; % V/s
% Ie = 1.5; % V/s
Ie = 2*ones(1, length(tvec));

[Iint, Vm, Theta] = MN_Calc(a, A, Ie, Tmax, nInt);

subplot(4, 5, 13)
yyaxis left; plot(tvec, Vm) 
ylim([1.1*min(Vm), max(Vm)+ovhd])
yyaxis right; plot(tvec, Ie) 
title('M: Tonic Bursting')

%% N - PHASIC BURSTING %%
Tmax = 0.5;
tvec = 0:dT:Tmax;
LT = length(tvec);
a = 5; % s^-1
A = [10, -0.6]; % V/s
% Ie = 1.5; % V/s
Ie = 1.5*ones(1, length(tvec));

[Iint, Vm, Theta] = MN_Calc(a, A, Ie, Tmax, nInt);

subplot(4, 5, 14)
yyaxis left; plot(tvec, Vm) 
ylim([1.1*min(Vm), max(Vm)+ovhd])
yyaxis right; plot(tvec, Ie) 
title('N: Phasic Bursting')

%% O - REBOUND BURST %% 
Tmax = 1;
tvec = 0:dT:Tmax;
LT = length(tvec);
a = 5; % s^-1
A = [10, -0.6]; % V/s
% Ie = 1.5; % V/s
Ie = ones(1, length(tvec));
Ie(1:LT/3) = 0;
Ie(LT/3+1:2*LT/3) = -3.5;
Ie(2*LT/3+1:end) = 0;

[Iint, Vm, Theta] = MN_Calc(a, A, Ie, Tmax, nInt);

subplot(4, 5, 15)
yyaxis left; plot(tvec, Vm) 
ylim([1.1*min(Vm), max(Vm)+ovhd])
yyaxis right; plot(tvec, Ie) 
ylim([-4.5, 1])
ylabel('external current (V/s)')
title('O: Rebound Burst')

%% P - MIXED MODE %% 
Tmax = 0.5;
tvec = 0:dT:Tmax;
LT = length(tvec);
a = 5; % s^-1
A = [5, -0.3]; % V/s
% Ie = 1.5; % V/s
Ie = 2*ones(1, length(tvec));

[Iint, Vm, Theta] = MN_Calc(a, A, Ie, Tmax, nInt);

subplot(4, 5, 16)
yyaxis left; plot(tvec, Vm) 
ylim([1.1*min(Vm), max(Vm)+ovhd])
ylabel('membrane potential (mV)')
yyaxis right; plot(tvec, Ie) 
title('P: Mixed Mode')
xlabel('time (s)');


%% Q - AFTERPOTENTIALS %% 
Tmax = 0.2;
tvec = 0:dT:Tmax;
LT = length(tvec);
a = 0; % s^-1
A = [0, 0]; % V/s
% Ie = 1.5; % V/s
Ie = ones(1, length(tvec));
Ie(1:LT/2) = 2;
Ie(LT/2+1:end) = 0;

[Iint, Vm, Theta] = MN_Calc(a, A, Ie, Tmax, nInt);

subplot(4, 5, 17)
yyaxis left; plot(tvec, Vm) 
ylim([1.1*min(Vm), max(Vm)+ovhd])
yyaxis right; plot(tvec, Ie) 
ylim([-1, 3])
title('Q: Afterpotentials')
xlabel('time (s)');

%% R - BASAL BISTABILITY %%
Tmax = 0.2;
tvec = 0:dT:Tmax;
LT = length(tvec);
a = 0; % s^-1
A = [8, -0.1]; % V/s
% Ie = 1.5; % V/s
Ie = ones(1, length(tvec));
Ie(1:LT/4) = 5;
Ie(LT/4+1:2*LT/4) = 0;
Ie(2*LT/4+1:3*LT/4) = 5;
Ie(3*LT/4+1:end) = 0;

[Iint, Vm, Theta] = MN_Calc(a, A, Ie, Tmax, nInt);

subplot(4, 5, 18)
yyaxis left; plot(tvec, Vm) 
ylim([1.1*min(Vm), max(Vm)+ovhd])
yyaxis right; plot(tvec, Ie) 
ylim([-1, 6])
title('R: Basal Bistability')
xlabel('time (s)');

%% S - PREFERRED FREQUENCY %% 
Tmax = 0.8;
tvec = 0:dT:Tmax;
LT = length(tvec);
a = 5; % s^-1
A = [-3, 0.5]; % V/s
% Ie = 1.5; % V/s
Ie = ones(1, length(tvec));
Ie(1:LT/8) = 5;
Ie(LT/8+1:2*LT/8) = 0;
Ie(2*LT/8+1:3*LT/8) = 4;
Ie(3*LT/8+1:4*LT/8) = 0;
Ie(4*LT/8+1:5*LT/8) = 5;
Ie(5*LT/8+1:6*LT/8) = 0;
Ie(6*LT/8+1:7*LT/8) = 4;
Ie(7*LT/8+1:end) = 0;

[Iint, Vm, Theta] = MN_Calc(a, A, Ie, Tmax, nInt);

subplot(4, 5, 19)
yyaxis left; plot(tvec, Vm) 
ylim([1.1*min(Vm), max(Vm)+ovhd])
yyaxis right; plot(tvec, Ie) 
ylim([-1, 6])
title('S: Preferred Frequency')
xlabel('time (s)');

%% T - SPIKE LATENCY %% 
Tmax = 0.05;
tvec = 0:dT:Tmax;
LT = length(tvec);
a = -80; % s^-1
A = [0, 0]; % V/s
% Ie = 1.5; % V/s
Ie = ones(1, length(tvec));
Ie(1:LT/2) = 8;
Ie(LT/2+1:end) = 0;

[Iint, Vm, Theta] = MN_Calc(a, A, Ie, Tmax, nInt);

subplot(4, 5, 20)
yyaxis left; plot(tvec, Vm) 
ylim([1.1*min(Vm), max(Vm)+ovhd])
yyaxis right; plot(tvec, Ie)
ylim([-1, 9])
title('T: Spike Latency')
xlabel('time (s)');
ylabel('external current (V/s)')
