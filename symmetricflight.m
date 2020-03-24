rho = 1.225;
S = 30;
c = 2.0569;
V = 90;
m = 6689.13;
g = 9.80665;
alpha0 = -0.0189;
alpha = 5.5*(pi/180) - alpha0;
W = m*g;
theta = alpha;
Z = -W * cos(theta);

C_z_alpha = -5.7434;
C_z_alphadot = -0.0035;
mu_c = m/(rho*S*c);
C_z_q = -5.6629;
C_m_alpha = -0.6615;
C_m_alphadot = 0.1780;
C_m_q = -8.7941;
K_Ysq = 1.3925;
C_x_u = -0.095;
C_z_0 = Z/(0.5*rho*V^2*S);
C_z_u = -0.3762;
C_m_u = 0.0699;
C_x_alpha = 0.4797;

%short period computations
Asp = 4 * (mu_c)^2 * (K_Ysq);
Bsp = -2 * (mu_c) * ((K_Ysq) * (C_z_alpha) + (C_m_alphadot) + (C_m_q));
Csp = (C_z_alpha) * (C_m_q) - 2 * (mu_c) * (C_m_alpha);
lambda_sp1 = ((- Bsp - 1j * sqrt(4* Asp * Csp - Bsp^2))/(2 * Asp))*(V/c);
lambda_sp2 = ((- Bsp + 1j * sqrt(4 * Asp * Csp - Bsp^2))/(2 * Asp))*(V/c);
w0_sp = (V/c) * sqrt(Csp/Asp);
dampratio_sp = Bsp/(2*sqrt(Asp*Csp))
Thalf_sp = - (0.693/real(lambda_sp1)) * (c/V);

%phugoid computations
Aph = 2 * (mu_c) * ((C_z_alpha) * (C_m_q) - 2 * (mu_c) * (C_m_alpha));
Bph = 2 * (mu_c) * ((C_x_u) * (C_m_alpha) - (C_m_u) * (C_x_alpha)) + ...
(C_m_q) * ((C_z_u) * (C_x_alpha) - (C_x_u) * (C_z_alpha));
Cph = (C_z_0) * ((C_m_u) * (C_z_alpha) - (C_z_u) * (C_m_alpha));
lambda_ph1 = ((- Bph - 1j * sqrt(4* Aph * Cph - Bph^2))/(2 * Aph))*(V/c);
lambda_ph2 = ((- Bph + 1j * sqrt(4* Aph * Cph - Bph^2))/(2 * Aph))*(V/c);
w0_ph = (V/c) * sqrt(Cph/Aph);
dampratio_ph = Bph/(2*sqrt(Aph*Cph))
Thalf_ph = - (0.693/real(lambda_ph1)) * (c/V);



