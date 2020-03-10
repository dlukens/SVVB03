%INPUT constants + Vc, hp, mf1, mf2, TAT 
Vc      = [251 221 190 161 134 121] / 1.9438444924574;          %[m/s]
hp      = [5030 5030 5020 5040 5040 5030] * 0.3048;             %[m]
mf1     = [770 638 533 444 420 420]/(3600*2.2046226218488);     %[kg/s]
mf2     = [806 670 579 488 450 450]/(3600*2.2046226218488);     %[kg/s]
TAT     = [10.2 7.8 5 4 2.2 1.5];                               %[K]
gamma   = 1.4;


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

for idx = 1:length(Vc)
    %calculating Mach number
    p = p0*(1+(lambda*hp(idx))/Temp0)^(g/(lambda*R));
    M(idx) = sqrt(2/(gamma-1)*((1+p0/p*((1+((gamma-1)*rho0*Vc(idx)^2)/(2*gamma*p0))^(gamma/(gamma-1))-1))^((gamma-1)/gamma)-1));
    %calculating temperature
    TISA = Temp0+lambda*hp(idx);
    Tm = TAT(idx) + TISA;
    T(idx) = Tm/(1+(gamma-1)/2*M(idx)^2);
    deltaT(idx) = T(idx)-TISA;
    %calculating Vt and Ve
    a = sqrt(gamma*R*T(idx));
    Vt(idx) = M(idx)*a;
    rho = p/(R*T(idx));
    Ve(idx) = Vt(idx)*sqrt(rho/rho0);
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

