%Run CG.m first!

%reference data
Vt1 = 89.6389656024862;
Vt2 = 89.7370024890751;

rho1 = 1.04229237795367;
rho2 = 1.03997427943693;

t1 = 51*60 + 2;
t2 = 52*60 + 46;

S = 30;

W1 = 9.80665*(AC_zerofuel(2) + F_init - FU(t1))/lb_unit;
W2 = 9.80665*(AC_zerofuel(2) + F_init - FU(t2))/lb_unit;

d_e1 = RefData.flightdata.delta_e.data(idx(t1));
d_e2 = RefData.flightdata.delta_e.data(idx(t2));

CN_1 = W1/(0.5*rho1*Vt1^2*S);
CN_2 = W2/(0.5*rho2*Vt2^2*S);

XCG_1 = (AC_zerofuel(3) + F_moment(AC_fuel(2) - FU(t1)))/(AC_zerofuel(2) + (F_init - FU(t1)));
XCG_2 = (AC_zerofuel(3) + F_moment(AC_fuel(2) - FU(t2)))/(AC_zerofuel(2) + (F_init - FU(t2)));

Cm_d = 1/(d_e2 - d_e1)*(CN_2*XCG_2 - CN_1*XCG_1)/MAC;


function cgx = F_moment(mass)
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