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

%% Stability of the Jacobian matrix (Question 3)
syms V w
dV_dt = (1/C)*(gCa*(0.5*(1+tanh((V-v1)/v2)))*(VCa-V) + gK*w*(VK-V) + gL*(VL-V)) + Iext;
dw_dt = phi*((0.5*(1+tanh((V-v3)/v4)))-w)*cosh((V-v3)/(2*v4));

JSymbolic = jacobian([dV_dt, dw_dt],[V,w]);
V = V_eq;
w = w_eq;
Jmatrix = zeros(2,2);
Jmatrix(1,1) = subs(JSymbolic(1,1));
Jmatrix(1,2) = subs(JSymbolic(1,2));
Jmatrix(2,1) = subs(JSymbolic(2,1));
Jmatrix(2,2) = subs(JSymbolic(2,2));

eigenValues = eig(Jmatrix);
fprintf('The eigen values are %d and %d \n', eigenValues(1), eigenValues(2));

%% Generating an action potential using MLE (Question 5)
options = odeset('RelTol',1e-3,'AbsTol',1e-6, 'refine',5, 'MaxStep', 1);

Iext = 0;
tSpan = [0, 300];
initial = [0,w_eq];

phi = 0.01;
[t1, S1] = ode15s(@(t,S)morris_lecar_ddt(t,S),tSpan, initial, options);

phi = 0.02;
[t2, S2] = ode15s(@(t,S)morris_lecar_ddt(t,S),tSpan, initial, options);

phi = 0.04;
[t3, S3] = ode15s(@(t,S)morris_lecar_ddt(t,S),tSpan, initial, options);

figure;
hold on;
plot(t1,S1(:,1));
plot(t2, S2(:,1));
plot(t3, S3(:,1));
xlabel('Time(in ms)');
ylabel('Volatage(in mV)');
title('Action potentials with different \phi');
legend('\phi = 0.01','\phi = 0.02','\phi = 0.04');
grid on;

figure;
hold on
Vnc = @(V) (Iext - gCa*(0.5*(1+tanh((V-v1)/v2)))*(V-VCa) - gL*(V-VL))/(gK*(V-VK));
wnc = @(V) (0.5*(1+tanh((V-v3)/v4)));
fplot(@(V) Vnc(V), [-80 100],'k');
fplot(@(V) wnc(V), [-80 100],'k');

plot(S1(:,1),S1(:,2));
plot(S2(:,1),S2(:,2));
plot(S3(:,1),S3(:,2));
xlabel('V(in mV)');
ylabel('w');
ylim([0,1]);
title('Phase Plane Plot(MLE)');
legend('\phi = 0.01','\phi = 0.02','\phi = 0.04');
grid on;

%% Depolarisation Threshold (Question 6)
Iext = 0;
tSpan = [0, 300];
phi=0.04;
initialV=linspace(-14.2,-13.8,401);  %to get threshold for phi = 0.04
%phi=0.01;
%initialV=linspace(-15.56,-15.48,401);    %to get threshold for phi = 0.01
max_V = zeros(401);

%phi = 0.01;
%[t1, S1] = ode15s(@(t,S)morris_lecar_ddt(t,S),tSpan, initial1, options);
%[t2, S2] = ode15s(@(t,S)morris_lecar_ddt(t,S),tSpan, initial2, options);


flag=1;
for n = 1:401 
    [t,S] = ode15s(@(t,S)morris_lecar_ddt(t,S),tSpan, [initialV(n),w_eq], options);
    max_V(n) = max(S(:,1));
    if max_V(n) >= 0 && flag==1
        fprintf("Threshold is (%f)",initialV(n))
        flag=0;
    end
end

% figure;
% hold on
% Vnc = @(V) (Iext - gCa*(0.5*(1+tanh((V-v1)/v2)))*(V-VCa) - gL*(V-VL))/(gK*(V-VK));
% wnc = @(V) (0.5*(1+tanh((V-v3)/v4)));
% fplot(@(V) Vnc(V), [-80 100],'k');
% fplot(@(V) wnc(V), [-80 100],'k');
% 
% plot(S1(:,1),S1(:,2));
% plot(S2(:,1),S2(:,2));

figure;
hold on
plot(initialV,max_V);

%% Response to higher Iext (Question 7)
Iext = 86;

figure;
hold on
Vnc1 = @(V) (Iext - gCa*(0.5*(1+tanh((V-v1)/v2)))*(V-VCa) - gL*(V-VL))/(gK*(V-VK));
wnc1 = @(V) (0.5*(1+tanh((V-v3)/v4)));
fplot(@(V) Vnc1(V), [-80 100]);
fplot(@(V) wnc1(V), [-80 100]);
xlabel('V(in mV)');
ylabel('w');
title('Phase Plane Plot(MLE)');

% Finding equilibrium points using MATLAB
syms V w
Vnc1_eqn = (1/C)*(Iext - gCa*(0.5*(1+tanh((V-v1)/v2)))*(V-VCa) - gK*w*(V-VK) - gL*(V-VL)) == 0;
wnc1_eqn = (0.5*(1+tanh((V-v3)/v4)) - w) == 0;
eq_pt_1 = solve([Vnc1_eqn, wnc1_eqn], [V, w]);

V_eq1 = double(eq_pt_1.V);
w_eq1 = double(eq_pt_1.w);

plot(V_eq1, w_eq1, 'k+', 'linewidth', 2);
text(V_eq1, w_eq1, ['(' num2str(round(V_eq1,5)) ',' num2str(round(w_eq1,5)) ')']);
grid on;

fprintf('The equilibrium point is located at (%d,%d) \n', V_eq1, w_eq1);

%------------------
tSpan = [0,300]
initial = [V_eq,w_eq];
[t1, S1] = ode15s(@(t,S)morris_lecar_ddt(t,S),tSpan, initial, options);
initial = [V_eq1,w_eq1];
[t2, S2] = ode15s(@(t,S)morris_lecar_ddt(t,S),tSpan, initial, options);
initial = [-27.9,0.17];
[t3, S3] = ode15s(@(t,S)morris_lecar_ddt(t,S),tSpan, initial, options);
% figure;
% hold on

plot(S1(:,1),S1(:,2));
plot(S2(:,1),S2(:,2));
plot(S3(:,1),S3(:,2));

end
function dS = HH_f_na(t,S)
    global C;
    global Iext;
    global gK;
    global gNa;
    global gL;
    global VK;
    global VNa;
    global VL;
    global eps;
    global f;
    
    V = S(1);
    n = S(2);
    m = S(3);
    h = S(4);
    
    alphan =  -0.01 * (V + eps + 50)/(exp(-(V + eps + 50)/10)-1);
    alpham =  -0.1 * (V + eps + 35)/(exp(-(V + eps + 35)/10)-1);
    alphah = 0.07 * exp(-(V + 60)/20);
    betan = 0.125 * exp(-(V + 60)/80);
    betam = 4 * exp(-(V + 60)/18);
    betah = 1/(exp(-(V + 30)/10) + 1);
    
    ddt_V = (1/C) * (Iext - gK * n^4 * (V - VK) - gNa * (1-f) * m^3 * h * (V-VNa) - gNa * f * m^3 * (V-VNa)  - gL * (V - VL));
    ddt_n = alphan * (1 - n) - betan * n;
    ddt_m = alpham * (1 - m) - betam * m;
    ddt_h = alphah * (1 - h) - betah * h;
    
    dS = [ddt_V; ddt_n; ddt_m; ddt_h];
end

function dS = HH_reduced(t,S)
    global C;
    global Iext;
    global gK;
    global gNa;
    global gL;
    global VK;
    global VNa;
    global VL;
    global eps;
    global f;
    global hconst;
    global Vr;
    V = S(1);
    n = S(2);
    
    alphan =  -0.01 * (V + eps + 50)/(exp(-(V + eps + 50)/10)-1);
    alpham =  -0.1 * (Vr + eps + 35)/(exp(-(Vr + eps + 35)/10)-1);
    alphah = 0.07 * exp(-(V + 60)/20);
    betan = 0.125 * exp(-(V + 60)/80);
    betam = 4 * exp(-(Vr + 60)/18);
    betah = 1/(exp(-(V + 30)/10) + 1);
    
    mInf = alpham/(alpham + betam);
    
    ddt_V = (1/C) * (Iext - gK * n^4 * (V - VK) - gNa * (1-f) * mInf^3 * hconst * (V-VNa) - gNa * f * mInf^3 * (V-VNa)  - gL * (V - VL));
    ddt_n = alphan * (1 - n) - betan * n;
    
    dS = [ddt_V; ddt_n];
end

function dS = HH_reduced_m(t,S,nInf,hInf)
    global C;
    global Iext;
    global gK;
    global gNa;
    global gL;
    global VK;
    global VNa;
    global VL;
    global eps;
    global Vr;
    
    V = S(1);
    m = S(2);
    
    alpham =  -0.1 * (V + eps + 35)/(exp(-(V + eps + 35)/10)-1);
    betam = 4 * exp(-(V + 60)/18);
    
    ddt_V = (1/C) * (Iext - gK * nInf^4 * (V - VK) - gNa * m^3 * hInf * (V - VNa) - gL * (V - VL));
    ddt_m = alpham * (1 - m) - betam * m;
    
    dS = [ddt_V; ddt_m];
end

function dS = hodgkin_huxley_ddt(t,S)
    global C;
    global Iext;
    global gK;
    global gNa;
    global gL;
    global VK;
    global VNa;
    global VL;
    global eps;
    
    V = S(1);
    n = S(2);
    m = S(3);
    h = S(4);
    
    alphan =  -0.01 * (V + eps + 50)/(exp(-(V + eps + 50)/10)-1);
    alpham =  -0.1 * (V + eps + 35)/(exp(-(V + eps + 35)/10)-1);
    alphah = 0.07 * exp(-(V + 60)/20);
    betan = 0.125 * exp(-(V + 60)/80);
    betam = 4 * exp(-(V + 60)/18);
    betah = 1/(exp(-(V + 30)/10) + 1);
    
    ddt_V = (1/C) * (Iext - gK * n^4 * (V - VK) - gNa * m^3 * h * (V - VNa) - gL * (V - VL));
    ddt_n = alphan * (1 - n) - betan * n;
    ddt_m = alpham * (1 - m) - betam * m;
    ddt_h = alphah * (1 - h) - betah * h;
    
    dS = [ddt_V; ddt_n; ddt_m; ddt_h];
end
%% Morris Lecar dynamics equation solver
function dS = morris_lecar_ddt(t,S)

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

%locally define state variables:
V=S(1);
w=S(2);

%local functions:
m_inf = (0.5*(1+tanh((V-v1)/v2)));
w_inf = (0.5*(1+tanh((V-v3)/v4)));

ddt_V = (1/C)*(gCa*m_inf*(VCa-V) + gK*w*(VK-V) + gL*(VL-V)+Iext);
ddt_w = phi*(w_inf-w)*cosh((V-v3)/(2*v4));

dS=[ddt_V; ddt_w];

end

