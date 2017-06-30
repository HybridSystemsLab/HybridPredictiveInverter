function xdot = f_inverter(x)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matlab M-file
% Project: Hybrid Predictive Inverter
%
% Name: f_inverter.m
%
% Description: Flow map
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global Vdc L R Cap omega simulateVDCnoise
% states
q   = x(1);
iL  = x(2);
vC  = x(3);
iR  = x(4);
vR  = x(5);
tau = x(6);

%% VDC noise 
if simulateVDCnoise
    Vdc = VDCnoise(tau);  
end

%% Flow map
qdot = 0;
iLdot = Vdc/L*q -R/L*iL - 1/L*vC;
vCdot = 1/Cap*(iL);
iRdot = -Cap*omega^2*vR;
vRdot = 1/Cap*iR;
taudot = 1;

xdot = [qdot; iLdot; vCdot; iRdot; vRdot; taudot];
end