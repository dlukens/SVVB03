clear;
clc;

MAC = 80.98; %[in]
X_MAC = 261.56; %[in]
lb_unit = 2.20462262185;
in_unit = 39.37007874;
F_init = 4100;
%F_init = 4050;

%Weights

Chipke1 = [131, 80];
Hans2 = [131, 102];
Gabriel3 = [214, 71];
Carlos4 = [214, 77];
Diego5 = [251, 76];
Lorenza6 = [251, 64];
Paula7  = [288, 74];
Luis8 = [288, 99];
Marta10 = [170, 60];

%%%Reference Data
% Chipke1 = [131, 95];
% Hans2 = [131, 92];
% Gabriel3 = [214, 66];
% Carlos4 = [214, 61];
% Diego5 = [251, 75];
% Lorenza6 = [251, 78];
% Paula7  = [288, 86];
% Luis8 = [288, 68];
% Marta10 = [170, 74];

%Matrix with weights[kg] and xcgdatum[in]
W = [Chipke1; Hans2; Gabriel3; Carlos4; Diego5; Lorenza6; Paula7; Luis8; Marta10];

%Convert weight[kg->lb]
W(:, 2) = W(:, 2) .* lb_unit;

%Add moment column[in-lb]
W(:, 3) = W(:, 1) .* W(:, 2);

%Payload (X_cgdatum[in], mass[lb], moment[in-lb]
AC_payload = sum(W,1);
AC_payload(1) = AC_payload(3)/AC_payload(2);

%AC balance from sheet (mass[lb], cg_datum[in], moment[in-lb])
AC_empty = [291.65, 9165, 2672953.5];

AC_zerofuel = AC_payload + AC_empty;
AC_zerofuel(1) = AC_zerofuel(3)/AC_zerofuel(2);

%Fuel load
AC_fuel = [285.5, F_init, 1170550];

%AC Ramp mass (X_cgdatum[in], mass[lb], moment[in-lb])
AC_ramp = AC_zerofuel + AC_fuel;
AC_ramp(1) = AC_ramp(3)/AC_ramp(2);

X_cg = AC_ramp(1) - X_MAC; %[in]
X_cg_chordpercent = X_cg/MAC * 100;
X_cg_m = X_cg/in_unit; %[m]

%AC mass[kg] at ramp
AC_rampmass_kg = AC_ramp(2)/lb_unit;

%%%%%%%%%%%Fuel as function of time%%%%%%%%%

FuelDat = load('FuelCG.dat');
FuelDat(:,2) = FuelDat(:,2) .* 100;
RefData = load('RefData.mat');
time = RefData.flightdata.time.data;
FU_left = RefData.flightdata.lh_engine_FU.data;
FU_right = RefData.flightdata.rh_engine_FU.data;

%Find function y = ax + b from data
global p;
p = polyfit(FuelDat(:,1), FuelDat(:,2), 1);

%Mass and cg of ac as function of time
t = 9;
AC_mass_t = F_init - FU(t);
CGdatum_t = (AC_zerofuel(3) + FU_moment(AC_fuel(2) - FU(t)))/(AC_zerofuel(2) + (F_init - FU(t)));

function cgx = FU_moment(mass)
    global p;
    cgx = p(1)*mass + p(2);
end

function i = idx(t)
    global time;
    i = find(time==round(t, 1));
end

function mass = FU(t)
    global FU_left;
    global FU_right;

    mass = FU_left(idx(t)) + FU_right(idx(t));
end

