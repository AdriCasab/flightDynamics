from math import *
import scipy.io as sio
import numpy as np

# Citation 550 - Linear simulation

# xcg = 0.25 * c

mat_contents = sio.loadmat('./data/matlab.mat')
struct = mat_contents['flightdata']
# print(struct.dtype) # gives an overview of the variables in the file
flightdata = struct[0, 0]

time_measurements = flightdata['time'][0][0][0][0]      # time in seconds

# Control surface deflections

d_a = np.radians(flightdata['delta_a'][0][0][0])        # aileron deflection [rad]
d_e = np.radians(flightdata['delta_e'][0][0][0])        # elevator deflection [rad]
d_r = np.radians(flightdata['delta_r'][0][0][0])        # rudder deflection [rad]
F_e = flightdata['column_fe'][0][0][0]                  # force on elevator control wheel [N]

# Aircraft response
roll_angle = np.radians(flightdata['Ahrs1_Roll'][0][0][0])
roll_rate = np.radians(flightdata['Ahrs1_bRollRate'][0][0][0])
yaw_rate = np.radians(flightdata['Ahrs1_bYawRate'][0][0][0])
pitch_rate = np.radians(flightdata['Ahrs1_bPitchRate'][0][0][0])


# Stationary flight condition, extract data

hp0 = flightdata['Dadc1_alt'][0][0][0] * 0.3048             # pressure altitude in the stationary flight condition [m]
V0 = flightdata['Dadc1_tas'][0][0][0] * 0.514444            # true airspeed in the stationary flight condition [m/sec]
V_computed = flightdata['Dadc1_cas'][0][0][0] * 0.514444    # Computer airspeed
alpha0 = np.radians(flightdata['vane_AOA'][0][0][0])        # angle of attack in the stationary flight condition [rad]
th0 = np.radians(flightdata['Ahrs1_Pitch'][0][0][0])        # pitch angle in the stationary flight condition [rad]
mach = flightdata['Dadc1_mach'][0][0][0]                    # mach number
Tstatic = flightdata['Dadc1_sat'][0][0][0] + 273.15         # static air temperature [K]
Ttotal = flightdata['Dadc1_tat'][0][0][0] + 273.15          # total air temperature [K]

'''
hp0 = 1527.  # pressure altitude in the stationary flight condition [m]
V0 = 129.9  # true airspeed in the stationary flight condition [m/sec] ffrom hp0
alpha0 = 0.8  # angle of attack in the stationary flight condition [rad]
th0 = .2  # pitch angle in the stationary flight condition [rad]
'''
# Aircraft mass

mass_pilots_coordinators = 80. + 102. + 101. # 95. + 92. + 74.
passenger1 = 86. # 68.
passenger2 = 76. # 78.
passenger3 = 81. # 75.
passenger4 = 76. # 66.
passenger5 = 84. # 61.
passenger6 = 85. # 86.
People = mass_pilots_coordinators + passenger1 + passenger2 + passenger3 + passenger4 + passenger5 + passenger6
initial_mass = 9165. / 2.20462 + 2700. / 2.20462 + People # 9165. / 2.20462 + 4050. / 2.20462 + People   # total mass [kg]

# fuel usage fuel flow doesn't work!!!!!!!

lh_engine_FU = flightdata['lh_engine_FU'][0][0][0] / 2.20462                # convert lbs to kg
rh_engine_FU = flightdata['rh_engine_FU'][0][0][0] / 2.20462                # convert lbs to kg
lh_engine_FMF = flightdata['lh_engine_FMF'][0][0][0] / 2.20462 / 3600       # convert lbs/hr to kg/s
rh_engine_FMF = flightdata['rh_engine_FMF'][0][0][0] / 2.20462 / 3600       # convert lbs/hr to kg/s


# aerodynamic properties
e = 0.8967405192324812 # 0.7447898545806619      # 0.8  # Oswald factor [ ]
CD0 = 0.02179637726372731 # 0.021153399475051692  # 0.04  # Zero lift drag coefficient [ ]
CLa = 4.761381912548548 # 4.611758223213205     # Slope of CL-alpha curve [ ]

# Longitudinal stability
Cma = -0.6302655337675444 # -0.7241557026400499   #-0.4969366132881471 #-0.70779  # longitudinal stabilty [ ]
Cmde = -1.316986174539227 # -1.5096732850136978  #-1.3653672995017823 #-1.4224  # elevator effectiveness [ ]

# Aircraft geometry

S = 30.00                   # wing area [m^2]
Sh = 0.2 * S                # stabiliser area [m^2]
Sh_S = Sh / S               # [ ]
lh = 0.71 * 5.968           # tail length [m]
c = 2.0569                  # mean aerodynamic cord [m]
lh_c = lh / c               # [ ]
b = 15.911                  # wing span [m]
bh = 5.791                  # stabilser span [m]
Ar = b ** 2 / S             # wing aspect ratio [ ]
Ah = bh ** 2 / Sh           # stabilser aspect ratio [ ]
Vh_V = 1                    # [ ]
ih = -2 * pi / 180          # stabiliser angle of incidence [rad]

# Constant values concerning atmosphere and gravity

rho0 = 1.2250               # air density at sea level [kg/m^3]
lambda_1 = -0.0065          # temperature gradient in ISA [K/m]
Temp0 = 288.15              # temperature at sea level in ISA [K]
R = 287.05                  # specific gas constant [m^2/sec^2K]
g = 9.81                    # [m/sec^2] (gravity constant)
pressure0 = 101325          # air pressure at sea level [Pa]
air_ratio = 1.4

# air density [kg/m^3]
'''
rho = rho0 * ((1 + (lambda_1 * hp0 / Temp0))) ** (-((g / (lambda_1 * R)) + 1))
W = m * g  # [N]       (aircraft weight)
'''
# Constant values concerning aircraft inertia
'''
muc = m / (rho * S * c)
mub = m / (rho * S * b)
'''
KX2 = 0.019
KZ2 = 0.042
KXZ = 0.002
KY2 = 1.25 * 1.114

# Aerodynamic constants

Cmac = 0  # Moment coefficient about the aerodynamic centre [ ]
CNwa = CLa  # Wing normal force slope [ ]
CNha = 2 * pi * Ah / (Ah + 2)  # Stabiliser normal force slope [ ]
depsda = 4 / (Ar + 2)  # Downwash gradient [ ]

# Lift and drag coefficient
'''
CL = 2 * W / (rho * V0 ** 2 * S)  # Lift coefficient [ ]
CD = CD0 + (CLa * alpha0) ** 2 / (pi * A * e)  # Drag coefficient [ ]
'''

# Stability derivatives

cmtc = -0.0064
#CX0 = W * sin(th0) / (0.5 * rho * V0 ** 2 * S)
CXu = -0.02792 # -0.095
CXa = +0.47966
CXadot = +0.08330
CXq = -0.28170
CXde = -0.03728

# CZ0 = -W * cos(th0) / (0.5 * rho * V0 ** 2 * S)
CZu = -0.37616
CZa = -5.74340
CZadot = -0.00350
CZq = -5.66290
CZde = -0.69612

Cmu = +0.06990
Cmadot = +0.17800
Cmq = -8.79415

CYb = -0.7500
CYbdot = 0
CYp = -0.0304
CYr = +0.8495
CYda = -0.0400
CYdr = +0.2300

Clb = -0.10260
Clp = -0.71085
Clr = +0.23760
Clda = -0.23088
Cldr = +0.03440

Cnb = +0.1348
Cnbdot = 0
Cnp = -0.0602
Cnr = -0.2061
Cnda = -0.0120
Cndr = -0.0939

# Dc = c / V0  # *d/dt
# Db = b / V0
# V = sqrt(W / (1 / 2. * CL * rho * S))


def get_time2idx(time):
    '''
    if it doesn't converge t = time - b, increase b
    :param time:
    :return: the correct index for the given input time
    '''
    t = int(time * 10. - 130.)
    idx = int(t * np.heaviside(t, 0))
    t0 = time_measurements[idx]
    t1 = time_measurements[idx+1]
    c = idx + 1
    while np.abs(time - t0) > np.abs(time - t1):
        c += 1
        t0 = t1
        t1 = time_measurements[c]

    return c - 1


def get_mass_update(time):

    t = get_time2idx(time)
    fuel_used = lh_engine_FU[t][0] + rh_engine_FU[t][0]  # sum of the mass flown from the 2 engines
    mass_update = initial_mass - fuel_used
    return mass_update


def get_ISA(time):
    '''
    simple function that returns ISA condition in gradient condition
    :param time: time from beginning of flight in seconds
    :return: T, rho, p in ISA condition
    '''
    t = get_time2idx(time)
    height = hp0[t][0]
    T = Temp0 + lambda_1 * height
    rho = rho0 * (T/Temp0) ** (-((g / (lambda_1 * R)) + 1))
    p = pressure0 * (T/Temp0) ** (-(g / (lambda_1 * R)))

    return T, rho, p


def get_stationary_condition(time):
    '''
    this functions returns actual values updated in time
    :param time: time in seconds, from start of flight test
    :return: V, mub, CL, muc, CX0, CZ0 values for symmetric e asymmetric matrices
    '''

    t = get_time2idx(time)          #in this case it works as index, *10 becasue the system works in 10 Hz

    # Stationary flight condition at a given point

    hp = hp0[t][0]                                  # pressure altitude in the stationary flight condition [m]
    V = V0[t][0]                                    # true airspeed in the stationary flight condition [m/sec]
    alpha = alpha0[t][0]                            # angle of attack in the stationary flight condition [rad] in stability axis system
    th = th0[t][0]                                  # pitch angle in the stationary flight condition [rad] in stability axis system
    mass_update = get_mass_update(time)
    T, rho, p = get_ISA(time)
    W = mass_update * g                             # [N] aircraft weight

    # varying values, concerning aircraft inertia
    '''try keep constant'''
    muc = mass_update / (rho * S * c)
    mub = mass_update / (rho * S * b)

    # Lift and drag coefficient

    CL = 2 * W / (rho * V ** 2 * S)                # Lift coefficient [ ]
    CD = CD0 + (CLa * alpha) ** 2 / (pi * Ar * e)  # Drag coefficient [ ]

    # Stabiblity derivatives

    CX0 = W * sin(th) / (0.5 * rho * V ** 2 * S)
    CZ0 = -W * cos(th) / (0.5 * rho * V ** 2 * S)

    return V, mub, CL, muc, CX0, CZ0


def get_arbitrary_variable(time, arr):
    """
    :param time: time in seconds
    :param arr:  array of values of interest read from .mat file
    :return: Value of variable of interest at time t
    """

    t = get_time2idx(time)
    return arr[t][0]


def get_mach(time):
    '''
    it uses ISA T, rho, p, the sea level T, rho, p, and V computed
    :param time: time from beginning of flight in seconds
    :return: mach number
    '''
    t = get_time2idx(time)
    Vc = V_computed[t][0]
    T, rho, p = get_ISA(time)
    member1 = 2. / (air_ratio - 1.)
    member2 = -1. + ((1. + ((air_ratio - 1.) * rho0 * (Vc**2.) / (2. * air_ratio * pressure0)))**(air_ratio / (air_ratio - 1.)))
    member3 = (1 + pressure0 / p * member2)**((air_ratio-1)/air_ratio) - 1
    M = sqrt(member1 * member3)
    return M


def get_Tstatic(time):

    Mach = get_mach(time)
    t = get_time2idx(time)
    Ttot = Ttotal[t][0]
    coeff = (air_ratio-1.)/2. * Mach**2.
    T_stat = Ttot/(1+coeff)
    return T_stat


def get_true_velocity(time):
    Mach = get_mach(time)
    T_stat = get_Tstatic(time)
    a = sqrt(air_ratio*R*T_stat)
    V_true = Mach*a
    return V_true


def get_equivalent_velocity(time):
    V_true = get_true_velocity(time)
    T, rho, p = get_ISA(time)
    V_eq = V_true * sqrt(rho/rho0)
    return V_eq


def get_reduced_equivalent_velocity(time):
    WS = 60500
    W = get_mass_update(time) * g
    ratio = WS / W
    V_eq = get_equivalent_velocity(time)
    v_reduced = V_eq * np.sqrt(ratio)
    return v_reduced


def get_reduced_Force(time):
    WS = 60500
    W = get_mass_update(time) * g
    ratio = WS / W
    Fe = get_arbitrary_variable(time, F_e)
    Fe_reduced = Fe * ratio
    return Fe_reduced
