%Cessna Citation II dimensions and coefficients
S= 30; %m2
cbar= 2.0569; %m
b= 15.911; %m
cd0= 0.0215;  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FROM FLIGHT DATA
cla= 4.4079; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FROM FLIGHT DATA
alpha0= -0.0189; % rad %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FROM FLIGHT DATA
e= 0.9521; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FROM FLIGHT DATA
m= 6689.13; %kg
g= 9.80665; %m/s^2
Var= 89; %m/s  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FROM PROCCESSED FLIGHT DATA
Vas= 92;  %m/s  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FROM PROCCESSED FLIGHT DATA
Vdr= 90; %m/s  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FROM PROCCESSED FLIGHT DATA
rho= 1.225;
kxx2= 0.019;
kyy2= 1.3925;
kzz2= 0.042;
kxz= 0.002;
mub= m/(rho*S*b);
%%%%%%%%% ASYMMETRIC STABILITY DERIVATIVES
%lateral force derivatives
cyb= -0.75;
cybdot= 0;
cyp= -0.0304;
cyr= 0.8495;
cyda= -0.04;
cydr= 0.23;
%roll moment derivatives
clb= -0.1026;
clp= -0.7108;
clr= 0.2376;
clda= -0.2309;
cldr= 0.0344;
%yaw moment derivatives
cnb= 0.1348;
cnbdot= 0;
cnp= -0.0602;
cnr= -0.2061;
cnda= -0.012;
cndr= -0.0939;

%%%%%%%%%%%% APERIODIC ROLL
alphaAR= (6.0*(pi/180)) -alpha0; %rad %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FROM PROCCESSED FLIGHT DATA
lambdaAR= (clp/(4*mub*kxx2))*(Var/b)

%%%%%%%%%%%% APERIODIC SPIRAL
cl= 2*m*g/(rho*Vas^2*S)
lambdaAS= ((2*cl*(clb*cnr-cnb*clr))/(clp*(cyb*cnr+4*mub*cnb)-cnp*(cyb*clr+4*mub*clb)))*(Vas/b)

%%%%%%%%%%%% DUTCH ROLL
alphaDR= (5.3*(pi/180)) -alpha0; %deg %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FROM PROCCESSED FLIGHT DATA
A_DR= -2*mub*kzz2
B_DR= 0.5*cnr
C_DR= -cnb
lambdaDR_1= ((- B_DR - 1j * sqrt(4* A_DR * C_DR - B_DR^2))/(2 * A_DR))*(Vdr/b)
lambdaDR_2= ((- B_DR + 1j * sqrt(4 * A_DR * C_DR - B_DR^2))/(2 * A_DR))*(Vdr/b)



