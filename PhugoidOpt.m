% Citation 550 - Linear simulation
% Aperiodic Roll
% xcg = 0.25*c

clear
clc

% Stationary flight condition
syms CZadotopt Cmadotopt CXuopt CXaopt CXqopt CZuopt CZaopt CZqopt Cmuopt Cmaopt Cmqopt CXdeopt CZdeopt Cmdeopt


load('FlightData.mat')

start = find(flightdata.time.data==3015); % For reference: 3225 For flight:3015
finish = find(flightdata.time.data==3135); % For reference: 3358 For flight: 3135

hp0    = 5030*0.3048;      	  % pressure altitude in the stationary flight condition [m]
V0     = flightdata.Dadc1_tas.data(start,1)*0.51444;            % true airspeed in the stationary flight condition [m/sec]
alpha0 = flightdata.vane_AOA.data(start,1)*pi/180 - (-0.0189);       	      % angle of attack in the stationary flight condition [rad]
th0    = flightdata.Ahrs1_Pitch.data(start,1)*pi/180;        % pitch angle in the stationary flight condition [rad]

% Aircraft mass
m      = 6720-(flightdata.lh_engine_FU.data(finish,1)+flightdata.rh_engine_FU.data(finish,1))*0.453592;         	  % mass [kg] Nuestro
%m      =  6689.13;         	  % mass [kg] Reference


% aerodynamic properties
e      = 0.9521;            % Oswald factor [ ] 
CD0    = 0.0215;            % Zero lift drag coefficient [ ]
CLa    = 4.4079;            % Slope of CL-alpha curve [ ]

% Longitudinal stability
Cma    = -0.6615;            % longitudinal stabilty [ ] Mahima: -0.7249 Ours: -0.6615
Cmde   = -1.4584;            % elevator effectiveness [ ] Mahima: -1.4968 Ours: -1.4584

% Aircraft geometry

S      = 30.00;	          % wing area [m^2]
Sh     = 0.2*S;           % stabiliser area [m^2]
Sh_S   = Sh/S;	          % [ ]
lh     = 0.71*5.968;      % tail length [m]
c      = 2.0569;	  % mean aerodynamic cord [m]
lh_c   = lh/c;	          % [ ]
b      = 15.911;	  % wing span [m]
bh     = 5.791;	          % stabilser span [m]
A_asym = b^2/S;           % wing aspect ratio [ ]
Ah     = bh^2/Sh;         % stabilser aspect ratio [ ]
Vh_V   = 1;		  % [ ]
ih     = -2*pi/180;       % stabiliser angle of incidence [rad]

% Constant values concerning atmosphere and gravity

rho0   = 1.2250;          % air density at sea level [kg/m^3] 
lambda = -0.0065;         % temperature gradient in ISA [K/m]
Temp0  = 288.15;          % temperature at sea level in ISA [K]
R      = 287.05;          % specific gas constant [m^2/sec^2K]
g      = 9.81;            % [m/sec^2] (gravity constant)

rho    = rho0*((1+(lambda*hp0/Temp0)))^(-((g/(lambda*R))+1));% [kg/m^3]  (air density)
rho = rho0;
W      = m*g;				                        % [N]       (aircraft weight)

% Constant values concerning aircraft inertia

muc    = m/(rho*S*c);
mub    = m/(rho*S*b);
KX2    = 0.019;
KZ2    = 0.042;
KXZ    = 0.002;
KY2    = 1.3925;

% Aerodynamic constants

Cmac   = 0;                     % Moment coefficient about the aerodynamic centre [ ]
CNwa   = CLa;   		        % Wing normal force slope [ ]
CNha   = 2*pi*Ah/(Ah+2);        % Stabiliser normal force slope [ ]
depsda = 4/(A_asym+2);               % Downwash gradient [ ]

% Lift and drag coefficient

CL = 2*W/(rho*V0*V0*S);               % Lift coefficient [ ]
CD = CD0 + (CLa*alpha0)^2/(pi*A_asym*e);  % Drag coefficient [ ]

% Stabiblity derivatives

CX0    = W*sin(th0)/(0.5*rho*V0^2*S); %x1
CXu    = -0.095  ;%x2
CXa    = 0.47966;%x3
CXadot = +0.08330;%x4
CXq    = -0.28170;%x5
CXde   = -0.03728;%x6

CZ0    = -W*cos(th0)/(0.5*rho*V0^2*S);%x7
CZu    = -0.37616;%x8
CZa    = -5.74340;%x9
CZadot = -0.00350;%x10
CZq    = -5.66290;%x11
CZde   = -0.69612;%x12

Cmu    = +0.06990;%x13
Cmadot = +0.17800;%x14
Cmq    = -8.79415;%x15

CYb    = -0.7500;%x16
CYbdot =  0     ;%x17
CYp    = -0.0304;%x18
CYr    = +0.8495;%x19
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

%% UNCORRELATED VARIABLES
SDList = [];
Cmalist = [];
init = -500;
last = 0;
step = 1;
i = 0;


% Cmde=	-1.0475;
% Cma=	-0.386;
% KY2=	1.00E-04;
% CZadot=	-476.7446;
% Cmadot=	10.1292;
% CXu	=-0.1152;
% CXa	=-1.3148;
% CZ0	=-0.5035;
% CXq	=-9.7885;
% CZu=	-0.538;
% CZa	=-9.2287;
% CX0	=-0.0493;
% CZq	=-16.23;
% Cmu	=0.088;
% Cmq	=-5.2958;
% CXde=	2.492;
% CZde=	-3.5501;


for idx = init:step:last
%    CZadot = idx
   i = i+1;

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
    
    
    C = eye(4);
    D = 0;
    
    t = flightdata.time.data(1,start:finish)-flightdata.time.data(1,start);
    
    sys_sym = ss(A_sym,B_sym,C,D);
    u_de = (flightdata.delta_e.data(start:finish,1)*pi/180)';
    
    y_sym = lsim(sys_sym,u_de,t);
    
    error1 = ((y_sym(:,1)+flightdata.Dadc1_tas.data(start,1)*0.51444)-flightdata.Dadc1_tas.data(start:finish,1)*0.51444)/(max(flightdata.Dadc1_tas.data(start:finish,1)*0.51444)-min(flightdata.Dadc1_tas.data(start:finish,1)*0.51444));
    % error1 = (y_sym(:,1)+flightdata.Dadc1_tas.data(start,1)*0.51444)-flightdata.Dadc1_tas.data(start:finish,1)*0.51444;
    variance1 = (dot(error1,error1))/length(y_sym(:,1));
    error2 = ((y_sym(:,2)+flightdata.vane_AOA.data(start,1)*pi/180)-flightdata.vane_AOA.data(start:finish,1)*pi/180)/(max(flightdata.vane_AOA.data(start:finish,1)*pi/180)-min(flightdata.vane_AOA.data(start:finish,1)*pi/180));
    % error2 = (y_sym(:,2)+flightdata.vane_AOA.data(start,1)*pi/180)-flightdata.vane_AOA.data(start:finish,1)*pi/180;
    variance2 = (dot(error2,error2))/length(y_sym(:,1));
%     variance2 = 0;
    error3 = ((y_sym(:,3)+flightdata.Ahrs1_Pitch.data(start,1)*pi/180)-flightdata.Ahrs1_Pitch.data(start:finish,1)*pi/180)/(max(flightdata.Ahrs1_Pitch.data(start:finish,1)*pi/180)-min(flightdata.Ahrs1_Pitch.data(start:finish,1)*pi/180));
    % error3 = (y_sym(:,3)+flightdata.Ahrs1_Pitch.data(start,1)*pi/180)-flightdata.Ahrs1_Pitch.data(start:finish,1)*pi/180;
    variance3 = (dot(error3,error3))/length(y_sym(:,1));
    error4 = ((y_sym(:,4)+flightdata.Ahrs1_bPitchRate.data(start,1)*pi/180)-flightdata.Ahrs1_bPitchRate.data(start:finish,1)*pi/180)/(max(flightdata.Ahrs1_bPitchRate.data(start:finish,1)*pi/180)-min(flightdata.Ahrs1_bPitchRate.data(start:finish,1)*pi/180));
    % error4 = (y_sym(:,4)+flightdata.Ahrs1_bPitchRate.data(start,1)*pi/180)-flightdata.Ahrs1_bPitchRate.data(start:finish,1)*pi/180;
    variance4 = (dot(error4,error4))/length(y_sym(:,1));
    
    SD = sqrt(variance1+variance2+variance3+variance4);
    SDList(i) = SD;
end

[SDmin, index] = min(SDList(:))
value = init+(index-1)*step

%% METHOD ONE
startcmde = -1;
endcmde = -0;
startcma = -20;
endcma = -15;
length1 = 10;
length2 = 10;
SDlist = [];
Cmdelist = [];
Cmalist = [];
SDMatrix=zeros(length1,length2);
i = 0;

for idx = startcmde:abs(startcmde-endcmde)/length1:endcmde
   CZu = idx; 
   i = i+1;
   j=0;
   Cmalist=[];
%%% Symmetric
for idx2 = startcma:abs(startcma-endcma)/length2:endcma
    j = j+1;
    SDlist = [];
    CZq = idx2;
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
    
    
    C = eye(4);
    D = 0;
    
    t = flightdata.time.data(1,start:finish)-flightdata.time.data(1,start);
    
    sys_sym = ss(A_sym,B_sym,C,D);
    u_de = (flightdata.delta_e.data(start:finish,1)*pi/180)';
    
    y_sym = lsim(sys_sym,u_de,t);
    
    error1 = ((y_sym(:,1)+flightdata.Dadc1_tas.data(start,1)*0.51444)-flightdata.Dadc1_tas.data(start:finish,1)*0.51444)/(max(flightdata.Dadc1_tas.data(start:finish,1)*0.51444)-min(flightdata.Dadc1_tas.data(start:finish,1)*0.51444));
    % error1 = (y_sym(:,1)+flightdata.Dadc1_tas.data(start,1)*0.51444)-flightdata.Dadc1_tas.data(start:finish,1)*0.51444;
    variance1 = (dot(error1,error1))/length(y_sym(:,1));
    error2 = ((y_sym(:,2)+flightdata.vane_AOA.data(start,1)*pi/180)-flightdata.vane_AOA.data(start:finish,1)*pi/180)/(max(flightdata.vane_AOA.data(start:finish,1)*pi/180)-min(flightdata.vane_AOA.data(start:finish,1)*pi/180));
    % error2 = (y_sym(:,2)+flightdata.vane_AOA.data(start,1)*pi/180)-flightdata.vane_AOA.data(start:finish,1)*pi/180;
    variance2 = (dot(error2,error2))/length(y_sym(:,1));
%     variance2 = 0;
    error3 = ((y_sym(:,3)+flightdata.Ahrs1_Pitch.data(start,1)*pi/180)-flightdata.Ahrs1_Pitch.data(start:finish,1)*pi/180)/(max(flightdata.Ahrs1_Pitch.data(start:finish,1)*pi/180)-min(flightdata.Ahrs1_Pitch.data(start:finish,1)*pi/180));
    % error3 = (y_sym(:,3)+flightdata.Ahrs1_Pitch.data(start,1)*pi/180)-flightdata.Ahrs1_Pitch.data(start:finish,1)*pi/180;
    variance3 = (dot(error3,error3))/length(y_sym(:,1));
    error4 = ((y_sym(:,4)+flightdata.Ahrs1_bPitchRate.data(start,1)*pi/180)-flightdata.Ahrs1_bPitchRate.data(start:finish,1)*pi/180)/(max(flightdata.Ahrs1_bPitchRate.data(start:finish,1)*pi/180)-min(flightdata.Ahrs1_bPitchRate.data(start:finish,1)*pi/180));
    % error4 = (y_sym(:,4)+flightdata.Ahrs1_bPitchRate.data(start,1)*pi/180)-flightdata.Ahrs1_bPitchRate.data(start:finish,1)*pi/180;
    variance4 = (dot(error4,error4))/length(y_sym(:,1));
    
    SD = sqrt(variance1+variance2+variance3+variance4);
    Cmalist = [Cmalist CZq];
    SDMatrix(i,j) = SD;
end
Cmdelist = [Cmdelist CZu];
idx
end

surf(Cmalist,Cmdelist,-SDMatrix);
[value, index] = min(SDMatrix(:));
[row, col] = ind2sub(size(SDMatrix), index);
value
Cmdelist(row)
Cmalist(col)

%% METHOD THREE

start1 = -1.11;
end1 = -1.09;
start2 = -0.01;
end2 = 0.01;
start3 = -0.6;
end3 = -0.4;
start4 = 0;
end4 = 5;
length1 = 10;
length2 = 10;
length3 = 10;
length4 = 1000;
SDlist = [];
C1list = [];
C2list = [];
C3list = [];
C4list = [];
i = 0;

for idx = start1:abs(start1-end1)/length1:end1
   CZu = idx; 
%%% Symmetric
    for idx2 = start2:abs(start2-end2)/length2:end2
        
        Cmu = idx2;
        
        
        for idx3 = start3:abs(start3-end3)/length3:end3
        k = 0; 
        CZ0 = idx3;
        
            for idx4 = start4:abs(start4-end4)/length4:end4
            
                CZq = idx4;
        
                k = k+1;
                CZ0 = idx3;
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


                C = eye(4);
                D = 0;

                t = flightdata.time.data(1,start:finish)-flightdata.time.data(1,start);

                sys_sym = ss(A_sym,B_sym,C,D);
                u_de = (flightdata.delta_e.data(start:finish,1)*pi/180)';

                y_sym = lsim(sys_sym,u_de,t);

                error1 = ((y_sym(:,1)+flightdata.Dadc1_tas.data(start,1)*0.51444)-flightdata.Dadc1_tas.data(start:finish,1)*0.51444)/(max(flightdata.Dadc1_tas.data(start:finish,1)*0.51444)-min(flightdata.Dadc1_tas.data(start:finish,1)*0.51444));
                % error1 = (y_sym(:,1)+flightdata.Dadc1_tas.data(start,1)*0.51444)-flightdata.Dadc1_tas.data(start:finish,1)*0.51444;
                variance1 = (dot(error1,error1))/length(y_sym(:,1));
                error2 = ((y_sym(:,2)+flightdata.vane_AOA.data(start,1)*pi/180)-flightdata.vane_AOA.data(start:finish,1)*pi/180)/(max(flightdata.vane_AOA.data(start:finish,1)*pi/180)-min(flightdata.vane_AOA.data(start:finish,1)*pi/180));
                % error2 = (y_sym(:,2)+flightdata.vane_AOA.data(start,1)*pi/180)-flightdata.vane_AOA.data(start:finish,1)*pi/180;
                variance2 = (dot(error2,error2))/length(y_sym(:,1));
            %     variance2 = 0;
                error3 = ((y_sym(:,3)+flightdata.Ahrs1_Pitch.data(start,1)*pi/180)-flightdata.Ahrs1_Pitch.data(start:finish,1)*pi/180)/(max(flightdata.Ahrs1_Pitch.data(start:finish,1)*pi/180)-min(flightdata.Ahrs1_Pitch.data(start:finish,1)*pi/180));
                % error3 = (y_sym(:,3)+flightdata.Ahrs1_Pitch.data(start,1)*pi/180)-flightdata.Ahrs1_Pitch.data(start:finish,1)*pi/180;
                variance3 = (dot(error3,error3))/length(y_sym(:,1));
                error4 = ((y_sym(:,4)+flightdata.Ahrs1_bPitchRate.data(start,1)*pi/180)-flightdata.Ahrs1_bPitchRate.data(start:finish,1)*pi/180)/(max(flightdata.Ahrs1_bPitchRate.data(start:finish,1)*pi/180)-min(flightdata.Ahrs1_bPitchRate.data(start:finish,1)*pi/180));
                % error4 = (y_sym(:,4)+flightdata.Ahrs1_bPitchRate.data(start,1)*pi/180)-flightdata.Ahrs1_bPitchRate.data(start:finish,1)*pi/180;
                variance4 = (dot(error4,error4))/length(y_sym(:,1));


                SD = sqrt(variance1+variance2+variance3+variance4);
                if k>1 && SD>SDlist(end)
                   break 
                end

                SDlist(end+1) = SD;
                if min(SDlist) == SD
                    C1list(end+1) = CZu;
                    C2list(end+1) = Cmu;
                    C3list(end+1) = CZ0;
                    C4list(end+1) = CZq;
                end
            
            end
        end
    end
idx
end

SD = min(SDlist)
CZu = C1list(end)
Cmu = C2list(end)
CZ0 = C3list(end)
CZq = C4list(end)


  