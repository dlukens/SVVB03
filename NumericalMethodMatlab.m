% Citation 550 - Linear simulation

% xcg = 0.25*c

% Stationary flight condition

hp0    = 5030*0.3048;      	  % pressure altitude in the stationary flight condition [m]
V0     = 90; %174.1*0.51444;            % true airspeed in the stationary flight condition [m/sec]
alpha0 = 5.5*pi/180;       	      % angle of attack in the stationary flight condition [rad]
th0    = 4.39444349167663*pi/180;        % pitch angle in the stationary flight condition [rad]

% Aircraft mass
m      = 6720;         	  % mass [kg] Nuestro
%m      =  6689.13;         	  % mass [kg] Reference


% aerodynamic properties
e      = 0.9521;            % Oswald factor [ ]
CD0    = 0.0215;            % Zero lift drag coefficient [ ]
CLa    = 4.4079;            % Slope of CL-alpha curve [ ]

% Longitudinal stability


% Aircraft geometry

S      = 30.00;	          % wing area [m^2]
Sh     = 0.2*S;           % stabiliser area [m^2]
Sh_S   = Sh/S;	          % [ ]5
lh     = 0.71*5.968;      % tail length [m]
c      = 2.0569;	  % mean aerodynamic cord [m]
lh_c   = lh/c;	          % [ ]
b      = 15.911;	  % wing span [m]
bh     = 5.791;	          % stabilser span [m]
A_asym      = b^2/S;           % wing aspect ratio [ ]
Ah     = bh^2/Sh;         % stabilser aspect ratio [ ]
Vh_V   = 1;		  % [ ]
ih     = -2*pi/180;       % stabiliser angle of incidence [rad]

% Constant values concerning atmosphere and gravity

rho0   = 1.2250;          % air density at sea level [kg/m^3] 
lambda = -0.0065;         % temperature gradient in ISA [K/m]
Temp0  = 288.15;          % temperature at sea level in ISA [K]
R      = 287.05;          % specific gas constant [m^2/sec^2K]
g      = 9.81;            % [m/sec^2] (gravity constant)

rho    = rho0*((1+(lambda*hp0/Temp0)))^(-((g/(lambda*R))+1));   % [kg/m^3]  (air density)
W      = m*g;				                        % [N]       (aircraft weight)

% Constant values concerning aircraft inertia

muc    = m/(rho*S*c);
mub    = m/(rho*S*b);
KX2    = 0.019;
KZ2    = 0.042;
KXZ    = 0.002;
KY2    = 1.25*1.114;

% Aerodynamic constants

Cmac   = 0;                     % Moment coefficient about the aerodynamic centre [ ]
CNwa   = CLa;   		        % Wing normal force slope [ ]
CNha   = 2*pi*Ah/(Ah+2);        % Stabiliser normal force slope [ ]
depsda = 4/(A_asym+2);               % Downwash gradient [ ]

% Lift and drag coefficient

CL = 2*W/(rho*V0^2*S);               % Lift coefficient [ ]
CD = CD0 + (CLa*alpha0)^2/(pi*A_asym*e);  % Drag coefficient [ ]

% Stabiblity derivatives

CX0    = W*sin(th0)/(0.5*rho*V0^2*S);
CXu    = -0.095  ;
CXa    = +0.47966;
CXadot = +0.08330;
CXq    = -0.28170;
CXde   = -0.03728;

CZ0    = -W*cos(th0)/(0.5*rho*V0^2*S);
CZu    = -0.37616;
CZa    = -5.74340;
CZadot = -0.00350;
CZq    = -5.66290;
CZde   = -0.69612;

Cmu    = +0.06990;
Cmadot = +0.17800;
Cmq    = -8.79415;

CYb    = -0.7500;
CYbdot =  0     ;
CYp    = -0.0304;
CYr    = +0.8495;
CYda   = -0.0400;
CYdr   = +0.2300;

Clb    = -0.10260;
Clp    = -0.71085;
Clr    = +0.23760;
Clda   = -0.23088;
Cldr   = +0.03440;

Cnb    =  +0.1348;
Cnbdot =   0     ;
Cnp    =  -0.0602;
Cnr    =  -0.2061;
Cnda   =  -0.0120;
Cndr   =  -0.0939;


%%% Symmetric

C1_sym = [-2*muc*c/(V0*V0), 0, 0, 0,;
      0, (CZadot - 2*muc)*c/V0, 0, 0;
      0, 0, -c/V0, 0;
      0, Cmadot*c/V0, 0, -2*muc*KY2*c*c/(V0*V0)];


C2_sym = [CXu/V0, CXa, CZ0, CXq*c/V0;
      CZu/V0, CZa, -CX0, (CZq + 2*muc)*c/V0;
      0,0,0,c/V0;
      Cmu/V0, Cma, 0, Cmq*c/V0];


C3_sym = [CXde;
      CZde;
      0;
      Cmde];

A_sym = -inv(C1_sym)*C2_sym;
B_sym = -inv(C1_sym)*C3_sym;

%%% Unsymmetric

C1_asym = [(CYbdot - 2*mub)*b/V0,0,0,0;
      0,-b/(2*V0),0,0;
      0,0,-2*mub*KX2*b*b/(V0*V0),2*mub*KXZ*b*b/(V0*V0);
      Cnbdot*b/V0,0,2*mub*KXZ*b*b/(V0*V0),-2*mub*KZ2*b*b/(V0*V0)];
  
  
C2_asym = [CYb, CL, CYp*b/(2*V0), (CYr - 4*mub)*b/(2*V0);
      0, 0, b/(2*V0), 0;
      Clb, 0, Clp*b/(2*V0), Clr*b/(2*V0);
      Cnb, 0, Cnp*b/(2*V0), Cnr*b/(2*V0)];
  
  
C3_asym = [CYda, CYdr;
       0,0;
       Clda, Cldr;
       Cnda, Cndr];
   
A_asym = -inv(C1_asym)*C2_asym;
B_asym = -inv(C1_asym)*C3_asym;

C = eye(4);
D = 0;

t = flightdata.time.data(1,30111:31741)-flightdata.time.data(1,30111);
sys_asym = ss(A_asym,B_asym,C,D);
sys_sym = ss(A_sym,B_sym,C,D);
y_asym = impulse(sys_asym,t);
%y_sym = impulse(sys_sym,t);
u_de = (flightdata.delta_e.data(30111:31741,1)*pi/180)';
%x0 = [V0,alpha0,th0,0];
y_sym = lsim(sys_sym,u_de,t);

%Assymetric Flight
% figure(1)
% subplot(4,1,1)
% plot(t,y_asym(:,1,1),t,y_asym(:,1,2))
% grid()
% 
% subplot(4,1,2)
% plot(t,y_asym(:,2,1),t,y_asym(:,2,2))
% grid()
% 
% subplot(4,1,3)
% plot(t,y_asym(:,3,1),t,y_asym(:,3,2))
% grid()
% 
% subplot(4,1,4)
% plot(t,y_asym(:,4,1),t,y_asym(:,4,2))
% grid()
% 
% %Symmetric Flight
% figure(2)
% subplot(4,1,1)
% plot(t,y_sym(:,1))
% grid()
% 
% subplot(4,1,2)
% plot(t,y_sym(:,2))
% grid()
% 
% subplot(4,1,3)
% plot(t,y_sym(:,3))
% grid()
% 
% subplot(4,1,4)
% plot(t,y_sym(:,4))
% grid()

%Validation

load('FlightData.mat')

figure(3)
%subplot(2,1,1)
%plot(t,y_sym(:,3)+alpha0,flightdata.time.data(1,30111:31741)-flightdata.time.data(1,30111),flightdata.Ahrs1_Pitch.data(30111:31741,1)*pi/180)
%plot(t,y_sym(:,1)+V0,flightdata.time.data(1,30111:31741)-flightdata.time.data(1,30111),flightdata.Dadc1_tas.data(30111:31741,1)*0.51444)

grid()

