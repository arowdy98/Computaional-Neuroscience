function morris_lecar

%declare model parameters
global C;
global gCa;
global VCa;
global gK;
global VK;
global gL;
global VL;
global v1;
global v2;
global v3;
global v4;
global phi;
global Iext;

%parameter values
C = 20 ; %microfarad/cm^2 
gCa=4.4; % millisiemens/ cm^2 
VCa=120; %millivolts
gK=8;% millisiemens/ cm^2 
VK=-84; %millivolts
gL=2;% millisiemens/ cm^2 
VL=-60;%millivolts
v1=-1.2; %millivolts
v2= 18 ; %millivolts
v3= 2 ; %millivolts
v4= 30; %millivolts
phi = 0.02; % per millisecond

Iext=0;

%% Generating nullclines, equilibrium points and trajectories (Question 2)
% generate nullclines
figure;
hold on
Vnc = @(V) (Iext - gCa*(0.5*(1+tanh((V-v1)/v2)))*(V-VCa) - gL*(V-VL))/(gK*(V-VK));
wnc = @(V) (0.5*(1+tanh((V-v3)/v4)));
fplot(@(V) Vnc(V), [-80 100]);
fplot(@(V) wnc(V), [-80 100]);
xlabel('V(in mV)');
ylabel('w');
title('Phase Plane Plot(MLE)');

% Finding equilibrium points using MATLAB
syms V w
Vnc_eqn = (1/C)*(Iext - gCa*(0.5*(1+tanh((V-v1)/v2)))*(V-VCa) - gK*w*(V-VK) - gL*(V-VL)) == 0;
wnc_eqn = (0.5*(1+tanh((V-v3)/v4)) - w) == 0;
eq_pt_0 = solve([Vnc_eqn, wnc_eqn], [V, w]);

V_eq = double(eq_pt_0.V);
w_eq = double(eq_pt_0.w);

plot(V_eq, w_eq, 'k+', 'linewidth', 2);
text(V_eq, w_eq, ['(' num2str(round(V_eq,5)) ',' num2str(round(w_eq,5)) ')']);
grid on;

fprintf('The equilibrium point is located at (%d,%d) \n', V_eq, w_eq);
% Simulate Trajectories
x = linspace(-80,100,100);
y = linspace(0,1,100);

[V_quiv,w_quiv] = meshgrid(x, y);

dV_dt = (1/C)*(gCa*(0.5*(1+tanh((V_quiv-v1)/v2))).*(VCa-V_quiv) + gK*w_quiv .*(VK-V_quiv) + gL*(VL-V_quiv)) + Iext;
dw_dt = phi*((0.5*(1+tanh((V_quiv-v3)/v4)))-w_quiv).*cosh((V_quiv-v3)/(2*v4));
quiver(x,y,dV_dt,dw_dt);
legend('V nullcline', 'w nullcline','Equilibrium Point','Trajectories');

end


