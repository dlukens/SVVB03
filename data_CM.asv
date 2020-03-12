%Run CG.m first!

t1 = 51*60 + 2;
t2 = 52*60 + 46;

W1 = 9.80665*(AC_zerofuel(2) + F_init - FU(t1));
W2 = 9.80665*(AC_zerofuel(2) + F_init - FU(t2));

d_e1 = RefData.flightdata.delta_e.data(idx(t1));
d_e2 = RefData.flightdata.delta_e.data(idx(t2));

%CN_1 = W1/(0.5*rho*Vt1^2*S);

%CN_2 = W2/(0.5*rho*Vt2^2*S);

%Cm_d = 1/(d_e2 - d_e1)*(CN_2*XCG_2 - CN_1*XCG_1)/MAC;




function i = idx(t)
    global time;
    i = find(time==round(t, 1));
end

function mass = FU(t)
    global FU_left;
    global FU_right;

    mass = FU_left(idx(t)) + FU_right(idx(t));
end