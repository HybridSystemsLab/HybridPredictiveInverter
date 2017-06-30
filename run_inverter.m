%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matlab M-file
% Project: Hybrid Predictive Inverter
%
% Filename: run_inverter
%
% Description: simulation of a Single-Phase DC/AC inverter controlled by 
%              a Hybrid Predictive Controller
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear
close all

%% Useful Constants
global omega R L Cap Vdc a b T rhoStar NumErr RL H Vout
global VoutMax deltaBar Levels lgc epsTilt GammaLim lambda

% Circuit parameters
R = 1;                  % resistance ||Ohms||
L = 2e-3;               % inductance ||Henry's||
Cap = 1.063e-3;         % capacitance ||Farads||
Vdc = 220;              % DC input voltage ||Volt||
f = 60;                 % frequency ||Hz||
omega = f*2*pi;         % angular freqency
RL = Inf;               % load resistance ||Ohms||
Levels = 3;             % inverter levels (see multilevel inverters)

% Tracking ellipse parameters
Vout = 100;             % desired output voltage ||Volt||
b = Vout;               % reference ellipse: vB axis
a = b*Cap*omega;        % reference ellipse: iB axis
ePerc = .1;             % tracking error percentage
rhoStar = (ePerc*Vout*Cap*omega)^2;
H = 1;                  % aspect ratio of the tracking ellipse
epsTilt = R*Cap/L;      % mixed term of function V
T = 0.001;             % prediction horizon ||s||

lgc = 1;                % used for the impulsive resets of z_r (see G.m)

% Controller constraints
lambda = .1;
k = abs(L*Cap*omega^2-1);
VoutMax = (Vdc/k-sqrt(rhoStar/((Cap*omega)^2 ...
    - (.5*epsTilt)^2)))*(k/(k + omega*R*Cap));
deltaBar = ((Cap*omega)^2-(.5*epsTilt)^2)*...
    ((Vdc-Vout*(omega*R*Cap + k))/k)^2;
GammaLim = (Vdc - omega*R*Cap*Vout)/k;

if Vout > VoutMax
    error('Vout <= Vout_Max: constraint not satisfied')
end
if deltaBar < rhoStar
    error('rhoStar <= deltaBar: constraint not satisfied')
end
if lambda <= 0 || lambda >= 1
    error('lambda \in (0,1): constraint not satisfied')
end

NumErr = 1e-10;

%% Initial conditions
global plotSolutions
plotSolutions = 0;  % Plot predicted trajectories at each jump
                    % (requires setting a breakpoint in g_inverter.m
if plotSolutions
    figure('Units', 'inches', 'Position',[10 10 12 10]), hold on;
end

global PredictionMethod
% FixedHorizon:     predicts solutions for all t \in [0, T]
% EventDetection:   predicts solutions until they hit the jump set
PredictionMethod = 'EventDetection';

global simulateVDCnoise simulateRefResets
simulateVDCnoise = 0;   % Adds two noise signals to the input voltage VDC 
simulateRefResets = 0;  % Simulates impulsive resets of the reference zr

q0 = 0;
thR = 0;
iR0 = omega*Cap*Vout*cos(thR);
vR0 = Vout*sin(thR);

aEll = sqrt(deltaBar)/sqrt(H);
bEll = sqrt(deltaBar)/(Cap*omega);
thetaEll = 1/2*atan(epsTilt/(H-(Cap*omega)^2));
th = 0;
iL0 = iR0 + sqrt(deltaBar/H)*cos(pi/2+1/100);
vC0 = vR0 + sqrt(deltaBar)/(Cap*omega)*sin(pi/2+1/100);

rho0 = rhoStar;
tau0 = 0;

x0 = [q0; iL0; vC0; iR0; vR0; tau0];

%% Simulation horizon
TSPAN = [0 .2];
JSPAN = [0 1e3];

% rule for jumps
% rule = 1 -> priority for jumps
% rule = 2 -> priority for flows
rule = 1;
solver = 'ode45';

%solver tolerances
options = odeset('RelTol',1e-5,'MaxStep',1e-5);
%% Simulation
[t, j, x] = HyEQsolver(@f_inverter, @g_inverter, @C_inverter, @D_inverter, ...
    x0,TSPAN,JSPAN,rule,options, solver);
%% Plot solution
postprocessing;