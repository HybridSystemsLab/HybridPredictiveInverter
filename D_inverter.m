function inside = D_inverter(x)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matlab M-file
% Project: Hybrid Predictive Inverter
%
% Name: D_inverter.m
%
% Description: Jump set
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global L R Cap Vdc omega rhoStar deltaBar H epsTilt lambda
% states
q   = x(1);
iL  = x(2);
vC  = x(3);
iR  = x(4);
vR  = x(5);
tau = x(6);

A = [0, -omega^2*Cap; 1/Cap, 0];
P = [H, epsTilt/2; epsTilt/2, (Cap*omega)^2];
ei = iL-iR;
ev = vC-vR;
e = [ei; ev];
ni = Vdc/L*q - R/L*iL + (L*Cap*omega^2-1)/L*vC;

V_B = e'*P*e;
dV_B = e'*(A'*P + P*A)*e + 2*e'*P*[ni;0];

if rhoStar <= V_B && V_B <= deltaBar
    inside = dV_B >= -lambda*R/L*V_B;
else
    inside = 0;
end
