file = 'RefStationaryData.xlsx';
% file = '20200305_V1.xlsx';

%INPUT constants + Vc, hp, mf1, mf2, TAT 

Vc      = readmatrix(file,'Range','E28:E33').' / 1.9438444924574;         %[m/s] CL-CD
hp      = readmatrix(file,'Range','D28:D33').' * 0.3048;                  %[m] CL-CD
mf1     = readmatrix(file,'Range','G28:G33').'/(3600*2.2046226218488);    %[kg/s] CL-CD
mf2     = readmatrix(file,'Range','H28:H33').'/(3600*2.2046226218488);    %[kg/s] CL-CD
TAT     = readmatrix(file,'Range','J28:J33').';                           %[K] CL-CD
mfused  = readmatrix(file,'Range','I28:I33').'/(2.2046226218488);         %[kg] CL-CD
alpha   = readmatrix(file,'Range','F28:F33').';                           %[deg] CL-CD

gamma   = 1.4;
dyn_visc= 1.7387652e-05;
Dinlet  = 0.686;                                                        %[m]
chord   = 2.0569;                                                       %[m]
S       = 30;                                                           %[m^2]
b       = 15.911;                                                       %[m^2]
A       = b^2/S;
g       = 9.81;                                                         %[m/sec^2] (gravity constant)

% total mass
minit   = 6689.2;                                  %[kg] reference
% minit   = 6719.9;                                   %[kg] test

% Constant values concerning atmosphere and gravity
rho0    = 1.2250;          % air density at sea level [kg/m^3] 
lambda  = -0.0065;         % temperature gradient in ISA [K/m]
Temp0   = 288.15;          % temperature at sea level in ISA [K]
R       = 287.05;          % specific gas constant [m^2/sec^2K]
g       = 9.80665;            % [m/sec^2] (gravity constant)
p0      = 101325;          % [Pa]

%OUTPUTS
M = zeros(1,length(Vc));
T = zeros(1,length(Vc));
deltaT = zeros(1,length(Vc));
Vt = zeros(1,length(Vc));
Ve = zeros(1,length(Vc));

%atmospheric conditions
p = zeros(1,length(Vc));
rho = zeros(1,length(Vc));

%% Tc

for idx = 1:length(Vc)
    %calculating Mach number
    p(idx) = p0*(1+(lambda*hp(idx))/Temp0)^(-g/(lambda*R));
    M(idx) = sqrt(2/(gamma-1)*((1+p0/p(idx)*((1+((gamma-1)*rho0*Vc(idx)^2)/(2*gamma*p0))^(gamma/(gamma-1))-1))^((gamma-1)/gamma)-1));
    %calculating temperature
    TISA = Temp0+lambda*hp(idx);
    Tm = TAT(idx) + 273.15;
    T(idx) = Tm/(1+((gamma-1)/2)*M(idx)^2);
    deltaT(idx) = T(idx)-TISA;
    %calculating Vt and Ve
    a = sqrt(gamma*R*T(idx));
    Vt(idx) = M(idx)*a;
    rho(idx) = p(idx)/(R*T(idx));
    Ve(idx) = Vt(idx)*sqrt(rho(idx)/rho0);
end

thrust_inputs = zeros(length(Vc),5);
thrust_inputs(:,1) = transpose(hp);
thrust_inputs(:,2) = transpose(M);
thrust_inputs(:,4) = transpose(mf1);
thrust_inputs(:,5) = transpose(mf2);
thrust_inputs(:,3) = transpose(deltaT);

%THRUST CALCULATION
%saving to matlab.dat
save matlab.dat thrust_inputs -ascii
%running thrust.exe
system('thrust.exe &');
%reading results
Thrust = importdata('thrust.dat');
%calculating thrust coefficients
Tc = zeros(1,length(Thrust));
for idx = 1:length(Thrust)
   Tc(idx) = (Thrust(idx,1)+Thrust(idx,2))/(0.5*rho(idx)*Vt(idx)^2*pi*((Dinlet)^2)/4);
end

%% CL AND CD

%drag coefficient
C_D = Tc*((pi*(Dinlet^2)/4)/S);
%lift coefficient
C_L = (minit*ones(1,length(Vc))-mfused)*g./(0.5*rho.*Vt.^2*S);

% plot(alpha,C_L)
% plot(C_L,C_D)
CL_vs_alpha = polyfit(alpha,C_L,1);
alpha0_deg = -CL_vs_alpha(2)/CL_vs_alpha(1);
alpha0_rad = alpha0_deg*pi/180;
CLa_deg = CL_vs_alpha(1);
CLa_rad = CL_vs_alpha(1)*180/pi
fitC_L = CLa_deg*alpha+CL_vs_alpha(2);

% plot(C_L.^2,C_D);
CD_vs_CLsq = polyfit(C_L.^2,C_D,1);
CD0 = CD_vs_CLsq(2)
e = 1/(pi*A*CD_vs_CLsq(1))
fitC_D = CD_vs_CLsq(1).*C_L.^2+CD_vs_CLsq(2);

Re = Vt*chord./dyn_visc;

%plots
figure(1)
clvsa = plot(alpha,C_L,'b');
clvsa.Marker='*';
hold on
fit1 = plot(alpha,fitC_L,'--');
dim = [0.5 0.12 0.5 0.11];
str = {'Mach range: 0.4148 - 0.2005','Re range: (1.6280 - 0.7848) e07'};
annotation('textbox',dim,'String',str,'FitBoxToText','on');
title('Lift Coefficient vs Angle of Attack');
xlabel('$$\alpha$$ [deg]','Interpreter','Latex');
ylabel('$$C_L$$ [ ]','Interpreter','Latex');
legend({'Flight Data','Linear Regression'},'Location','northwest')
grid on

figure(2)
clvscd = plot(C_L,C_D,'r');
clvscd.Marker='*';
hold on
fit2 = plot(C_L,fitC_D,'--');
dim = [0.5 0.12 0.5 0.11];
str = {'Mach range: 0.4148 - 0.2005','Re range: (1.6280 - 0.7848) e07'};
annotation('textbox',dim,'String',str,'FitBoxToText','on');
title('Drag Coefficient vs Lift Coefficient');
xlabel('$$C_L$$ [ ]','Interpreter','Latex');
ylabel('$$C_D$$ [ ]','Interpreter','Latex');
legend({'Flight Data','Lift Squared Regression'},'Location','northwest')
grid on



