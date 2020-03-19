import numpy as np
import scipy.io
from scipy import interpolate

FlightData = scipy.io.loadmat('FlightData.mat')
RefData = scipy.io.loadmat('RefData.mat')


MAC = 80.98 #[in]
X_MAC = 261.56 #[in]
lb2kg = 2.20462262185
m2in = 39.37007874

#######Payload Weights and location#####

F_init = 4100/lb2kg;
W = np.array([[131., 80.],
            [131., 102.],
            [214., 71.],
            [214., 77.],
            [251., 76.],
            [251., 64.],
            [288., 74.],
            [288., 99.],
            [170., 60.]])

W[:, 0] = W[:, 0]/m2in
Wb = np.array(W)
Wb[7][0] = 134./m2in

###Reference Data
# F_init = 4050/lb2kg;
# W = np.array([[131., 95.],
#             [131., 92.],
#             [214., 66.],
#             [214., 61.],
#             [251., 75.],
#             [251., 78.],
#             [288., 86.],
#             [288., 68.],
#             [170., 74.]])

# W[:, 0] = W[:, 0]/m2in
# Wb = np.array(W)
# Wb[7][0] = 134./m2in

#Add moment column
W = np.hstack((W, np.atleast_2d(W[:, 0]*W[:, 1]).T))
Wb = np.hstack((Wb, np.atleast_2d(Wb[:, 0]*Wb[:, 1]).T))

###FUNCTIONS####
def F_moment(mass):
    FuelCG = np.loadtxt('FuelCG.dat', delimiter=',')
    FuelCG[:,0] *= 1/lb2kg
    FuelCG[:,1] *= 100/m2in/lb2kg
    f = interpolate.interp1d(FuelCG[:,0], FuelCG[:,1])
    
    return f(mass)


#AIRCRAFT WEIGHTS####
#AC payload weight and moment. [kgf], [kgf m]
AC_payload = np.array((np.sum(W[:, 1]), np.sum(W[:, 2])))
AC_payloadb = np.array((np.sum(Wb[:, 1]), np.sum(Wb[:, 2])))

AC_empty = np.array((9165/lb2kg, 2672953.5/m2in/lb2kg))

AC_zerofuel = np.array(AC_empty + AC_payload)
AC_zerofuelb = np.array(AC_empty + AC_payloadb)

AC_fuel = np.array((F_init, F_moment(F_init)))

AC_ramp = AC_zerofuel + AC_fuel
AC_rampb = AC_zerofuelb + AC_fuel

###Fuel as funtion of time

F_usedt1 = 881/lb2kg
F_usedt2 = 910/lb2kg

x_CG1 = (AC_ramp[1] - F_moment(F_usedt1))/(AC_ramp[0] - F_usedt1)
x_CG2 = (AC_rampb[1] - F_moment(F_usedt2))/(AC_rampb[0] - F_usedt2)

diff = x_CG1 - x_CG2
