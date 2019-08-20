function [Iint, Vm, Theta, f] = MN_Calc(a, A, Ie, Tmax, nInt)

%%% Input Arguments: a, A, Ie, maxIter
%%%
%%%          a: simulation parameter, [a], typically 0 or 5 (1/s)
%%%          A: simulation parameter, [A1/C, A2/C] (V/s)
%%%          Ie: constant external current, [Ie/C], between 0-8 (V/s)
%%%          maxIter: max time, [Tmax], typically between 0.2-0.5 (s)
%%%          nInt: number of internal currents, [nInt] (unit)
%%%
%%% Output Arguments: Iint, Vm, Theta
%%%
%%%          Iint: internal currents, [Iint/C] (V/s)
%%%          Vm: membrane voltage, [Vm] (V)
%%%          Theta: instantaneous threshold, [Vm] (V)
%%%
%%%
%%% Written by:    Khaled Aboumerhi
%%% Date:          08/19/2019
%%% Institution:   Johns Hopkins University
%%% Contact Info:  aboumerhi@jhu.edu


b = 10; % s^-1
EL = -0.07; % V
Vr = -0.07; % resting membrane voltage (V)
G = 50; % s^-1
rtheta = -0.06; % V
itheta = -0.05; % V
R = [0, 1];
kval = [200, 20]; % s^-1

dT = 1E-4; % time step (s)
tvec = 0:dT:Tmax; 

Iint = zeros(nInt, length(tvec)); % internal currents (V/s)
dIint_dT = zeros(nInt, 1); 
Vm = zeros(1, length(tvec)); % membrane voltage (V)
Vm(1) = Vr;
Theta = zeros(1, length(tvec)); % instantaneous threshold
Theta(1) = itheta;

dVm_dT = 0;
dTheta_dT = 0;

for iter = 1:length(tvec)-1
    
    % time evolution of the state variables    
    for ii = 1:length(kval)
        dIint_dT(ii) = -kval(ii)*Iint(ii, iter);
    end
    dVm_dT = Ie(iter) + sum(Iint(:, iter)) - G*(Vm(iter) - EL);
    dTheta_dT = a*(Vm(iter) - EL) - b*(Theta(iter) - itheta);
    
    % update differential variables   
    dInt = dIint_dT*dT;
    dVm = dVm_dT*dT;
    dTheta = dTheta_dT*dT;
    
    Iint(:, iter + 1) = Iint(:, iter)  + dInt;
    Vm(iter + 1) = Vm(iter) + dVm;
    Theta(iter + 1) = Theta(iter) + dTheta;
    
    % update rule only applies if Vm < instantaneous threshold
    if Vm(iter + 1) > Theta(iter + 1)
        Iint(:, iter + 1) = R*Iint(:, iter + 1) + A;
        Vm(iter + 1) = Vr;
        Theta(iter + 1) = max(rtheta, Theta(iter+ 1));
    end

end

% Vm(end) = [];
Vm = Vm*10^(3); % V -> mV

% keyboard
% for t = 2:maxIter
%     
%     for n = 1:length(k)
%         Iint(t, n) = -k(n)*Iint(t-1, n);
%     end
%     
%     Vm(t) = Ie + sum(Iint(t-1, :)) - G*(Vm(t-1) - EL);
%     Theta(t) = a*(Vm(t-1) - EL) - b*(Theta(t-1) - thetaMax);
%     
%     if Vm(t) == Theta(t)
%         
%         for n = 1:length(k)
%             Iint(t, n) = R(n)*Iint(t, n) + A(n);
%         end
%         
%         Vm(t) = Vr;
%         Theta(t) = max(Thetar, Theta(t));
%         
%     end
% end
% 
% end
% 
