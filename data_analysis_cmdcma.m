% file = 'RefStationaryData.xlsx';
file = '20200305_V1.xlsx';

%INPUT constants + Vc, hp, mf1, mf2, TAT 

Vc      = readmatrix(file,'Range','E75:E76').' / 1.9438444924574;         %[m/s] Trim
hp      = readmatrix(file,'Range','D75:D76').' * 0.3048;                  %[m] Trim
mf1     = readmatrix(file,'Range','J75:J76').'/(3600*2.2046226218488);    %[kg/s] Trim
mf2     = readmatrix(file,'Range','K75:K76').'/(3600*2.2046226218488);    %[kg/s] Trim
TAT     = readmatrix(file,'Range','M75:M76').';                           %[K] Trim
mfused  = readmatrix(file,'Range','L75:L76').'/(2.2046226218488);         %[kg] Trim
alpha   = readmatrix(file,'Range','F75:F76').';                           %[deg] Trim
de_meas = readmatrix(file,'Range','G75:G76').';                           %[deg] Trim
Fe      = readmatrix(file,'Range','I75:I76').';                           %[deg] Trim

gamma   = 1.4;
Dinlet  = 0.691;                                                          %[m]
chord   = 2.0569;                                                         %[m]
S       = 30;                                                             %[m^2]
b       = 15.911;                                                         %[m^2]
A       = b^2/S;
g       = 9.81;                                                         %[m/sec^2] (gravity constant)

% total mass
% minit   = 6689.22;                                  %[kg] reference
minit   = 6719.9027;                                   %[kg] flight test

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

%standard values
mfs = 0.048; %kg/s per engine
ms = 60500/g; %kg

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
   Tc(idx) = (Thrust(idx,1)+Thrust(idx,2))/(0.5*rho(idx)*Vt(idx)^2*((Dinlet)^2));
end

%% CMDE
% xcg1 = 7.1447; %reference data
xcg1 = 7.1511; %flight data
% xcg2 = 7.1022; %reference data se mueve el delgado
% xcg2 = 7.0910; %reference data se mueve el gordo
% xcg2 = 7.0898; %flight data to 134
xcg2 = 7.0883; %flight data to 131

CN2 = ((minit-mfused(2))*g)/(0.5*rho(2)*Vt(2)^2*S);
CN1 = ((minit-mfused(1))*g)/(0.5*rho(1)*Vt(1)^2*S);

Cmde = -1/(de_meas(2)*pi/180-de_meas(1)*pi/180)*CN2*(xcg2-xcg1)/(chord)

%% CMA
alpha_der   = readmatrix(file,'Range','F59:F65').';                           %[deg] Trim
de_meas_der = readmatrix(file,'Range','G59:G65').';                           %[deg] Trim
slope_for_cma = polyfit(alpha_der,de_meas_der,1);

Cma = -slope_for_cma(1)*Cmde