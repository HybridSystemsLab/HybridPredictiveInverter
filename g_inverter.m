function xplus = g_inverter(x)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matlab M-file
% Project: Hybrid Predictive Inverter
%
% Name: g_inverter.m
%
% Description: Jump map
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global plotSolutions PredictionMethod simulateRefResets
global T omega Cap R L Vdc Levels lgc Vout
% states
q = x(1);
iL = x(2);
vC = x(3);
iR = x(4);
vR = x(5);
tau = x(6);
% simulation horizon
TPspan = [0 T];
% admissible values of q
q0 = linspace(-1,1,Levels);

%% Jump map: qplus
Tif = zeros(size(q0)); % Initialize time to impact function

% Compute qhat according to the map Ghat
qbar = (R*iR - (L*Cap*omega^2-1)*vC)/Vdc;
if iL-iR < -R*Cap/(2*L)*(vC-vR)
    qhat = q0(q0 >= qbar);
elseif iL-iR > -R*Cap/(2*L)*(vC-vR)
    qhat = q0(q0 <= qbar);
else
    qhat = q0;
end
% Compute qplus
switch PredictionMethod
    case 'FixedHorizon'
        options = odeset('RelTol',1e-5,'MaxStep',1e-4,'Events', []);
        if plotSolutions
            qhat = q0;
        end
        for i = 1:length(qhat)
            x0 = [qhat(i), iL, vC, iR, vR, 0];
            [tP, xP] = ode45(@(t, xP) f_inverter(xP), TPspan, x0, options);
            Tif(i) = 1;
            % Increment Iq if xP \in C
            while Tif(i) <= size(xP,1) && ~D_inverter(xP(Tif(i),:))
                Tif(i) = Tif(i) + 1;
            end
            if plotSolutions
                len(i) = Tif(i);
                solPlot{i}.t = tP;
                solPlot{i}.x = xP;
            end
        end
        qplus = qhat(Tif == max(Tif));
        
    case 'EventDetection'
        options = odeset('RelTol',1e-5,'MaxStep',1e-4,'Events', ...
            @odeJumpEvent);
        if plotSolutions
            qhat = q0;
        end
        for i = 1:length(qhat)
            x0 = [qhat(i), iL, vC, iR, vR, 0];
            if ~D_inverter(x0) 
                [tP, xP]=ode45(@(t, xP) f_inverter(xP),TPspan,x0,options);
                Tif(i) = tP(end);
            else
                tP = 0;
                xP = x0;
            end
            if plotSolutions
                len(i) = length(tP);
                solPlot{i}.t = tP;
                solPlot{i}.x = xP;
            end
        end
        qplus = qhat(Tif == max(Tif));
    otherwise
        qplus = q;
end

% If the jump map is set-valued, randomly select one qplus among G(x)
if length(qplus) > 1
    qplus = qplus(ceil(rand*length(qplus)));
end

% Plot predicted trajectories
if plotSolutions
    plotPredTrajectories;
    % (set breakpoint here)
end

% Simulate impulsive reference resets
if simulateRefResets
    if tau >= 0.05
        Vout = Vout - lgc*Vout/2;
        lgc = -lgc;
        theta = rand*2*pi;
        iRplus = Cap*omega*Vout*cos(theta);
        vRplus = Vout*sin(theta);
        tauplus = 0;
    else
        tauplus = tau;
        iRplus = iR;
        vRplus = vR;
    end
else
    iLplus = iL;
    vCplus = vC;
    iRplus = iR;
    vRplus = vR;
    tauplus = tau;
end

xplus = [qplus; iLplus; vCplus; iRplus; vRplus; tauplus];
end