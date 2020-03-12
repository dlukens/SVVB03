file = 'RefStationaryData.xlsx';

%INPUT constants + Vc, hp, mf1, mf2, TAT 
Vc      = readmatrix(file,'Range','E28:E33') / 1.9438444924574;          %[m/s]
hp      = readmatrix(file,'Range','D28:D33') * 0.3048;             %[m]
mf1     = readmatrix(file,'Range','G28:G33')/(3600*2.2046226218488);     %[kg/s]
mf2     = readmatrix(file,'Range','H28:H33')/(3600*2.2046226218488);     %[kg/s]
TAT     = readmatrix(file,'Range','J28:J33');                               %[K]
gamma   = 1.4;
Dinlet  = 0.686;                                                %[m]


% Constant values concerning atmosphere and gravity
rho0    = 1.2250;          % air density at sea level [kg/m^3] 
lambda  = -0.0065;         % temperature gradient in ISA [K/m]
Temp0   = 288.15;          % temperature at sea level in ISA [K]
R       = 287.05;          % specific gas constant [m^2/sec^2K]
g       = 9.81;            % [m/sec^2] (gravity constant)
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

for idx = 1:length(Vc)
    %calculating Mach number
    p(idx) = p0*(1+(lambda*hp(idx))/Temp0)^(g/(lambda*R));
    M(idx) = sqrt(2/(gamma-1)*((1+p0/p(idx)*((1+((gamma-1)*rho0*Vc(idx)^2)/(2*gamma*p0))^(gamma/(gamma-1))-1))^((gamma-1)/gamma)-1));
    %calculating temperature
    TISA = Temp0+lambda*hp(idx);
    Tm = TAT(idx) + TISA;
    T(idx) = Tm/(1+(gamma-1)/2*M(idx)^2);
    deltaT(idx) = T(idx)-TISA;
    %calculating Vt and Ve
    a = sqrt(gamma*R*T(idx));
    Vt(idx) = M(idx)*a;
    rho(idx) = p(idx)/(R*T(idx));
    Ve(idx) = Vt(idx)*sqrt(rho(idx)/rho0);
end

thrust_inputs = zeros(6,5);
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
   Tc(idx) = (Thrust(idx,1)+Thrust(idx,2))/(0.5*rho(idx)*Vt(idx)^2*Dinlet^2);
end

