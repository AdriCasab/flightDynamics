"""
Group 1203F1 Flight Dynamics Assignment
Symmetric Numerical Model - Characterizing and Simulating Eigenmotions
"""
from Cit_par import *
import matplotlib.pyplot as plt
import numpy as np
from numpy.linalg import inv
from eigen_methods import characterize_eigenmotions
import control


'Create State Space system'

def get_symmetric_state_space(time):
    '''
    from given initial time return state space for symmetric motion
    :param time: insert initial time in seconds at start of eigenmotion
    :return: matrix A and state space for symmetric motion
    '''
    V, mub, CL, muc, CX0, CZ0 = get_stationary_condition(time)

    'Derivation for state space'

    P = np.array([[-2 * muc * c / (V ** 2.), 0, 0, 0],
                  [0, (CZadot - 2 * muc) * c / V, 0, 0],
                  [0, 0, -c / V, 0],
                  [0, Cmadot * c / V, 0, -2 * muc * KY2 * (c / V) ** 2]])

    Q = np.array([[-CXu / V, -CXa, -CZ0, -CXq * (c / V)],
                  [-CZu / V, -CZa, CX0, -(CZq + 2 * muc) * (c / V)],
                  [0, 0, 0, -(c / V)],
                  [-Cmu / V, -Cma, 0, -Cmq * (c / V)]])

    R = np.array([[-CXde], [-CZde], [0], [-Cmde]])

    'Report derivation constants'

    C1 = P
    C2 = -Q
    C3 = -R

    'State space representation'

    A = inv(P) @ Q  # state matrix
    B = inv(P) @ R  # input matrix
    C = np.identity(A.shape[0])  # output matrix
    D = np.zeros((4, 1))  # direct matrix
    print(A)
    sys = control.ss(A, B, C, D)

    return A, sys


'Get data from flightest file'


def get_validation(t_sim):
    '''
        from given time return specified params
        :param t_sim: simulation time vector
        :return extract vector parameters (velocity, alpha, pitch, pitch_rate, elevator deflection) for t_sim
        '''
    vel = []
    alpha = []
    theta = []
    q = []
    delta_e = []

    for time in t_sim:
        vel.append(get_arbitrary_variable(time, V0))
        alpha.append(get_arbitrary_variable(time, alpha0))
        theta.append(get_arbitrary_variable(time, th0))
        q.append(get_arbitrary_variable(time, pitch_rate))
        delta_e.append(get_arbitrary_variable(time, d_e))
    return vel, alpha, theta, q, delta_e


'Simulation at given time for input taken from flightest data'


def response_system(t, sys, string, inputvec=0):
    if string == "step":
        y = control.step_response(sys, t, U=np.transpose(inputvec))[1]
    elif string == "forced":
        y = control.forced_response(sys, t , U=np.transpose(inputvec), X0=0.)[1]
    elif string == "impulse":
        y = control.impulse_response(sys, t, U=np.transpose(inputvec), X0=0.0)[1]
    elif string == "inputvec":
        y = control.matlab.lsim(sys, inputvec, t)[0]
        y = np.transpose(y)
    else:
        raise ValueError("Invalid input type specified")
    return y  # u, alpha, theta, q


def sim_sym_resp(emotion, ti, tf, t_plot=None, dt=0.1, plot=True, n_periods=1, d_e=None, delay=0.0):
    """
        :param emotion  : string that define what emotion is simulated
        :param ti       : time in seconds when eigenmotion starts
        :param dt       : time step in seconds for the simulation
        :param tf       : time/duration in seconds of the width input pulse
        :param plot     : boolean to show the output parameters
        :param n_periods: integer to simulate given number of periods if tf not specified
        :param d_e      : average deflection of elevator to simulate the input vector to feed to response_system()
        :param delay    : delay of the simulated input (real response takes more time) / tweaking
        :return time vector and output variables vector
        """
    A, sys = get_symmetric_state_space(ti)  # insert time in seconds
    print(characterize_eigenmotions(A))
    period = characterize_eigenmotions(A)['Period']
    t_short = min(abs(period))  # T_short_period[s]
    t_long = max(abs(period))   # T_phugoid[s]
    # move this line to 120 to make it step instead of pulse (here is more correct wrt to real delfection)
    deltas = np.append(np.zeros(int(delay/dt)), np.ones(int((tf-delay)/dt)) * d_e)

    # if tf is None:
    #     if emotion == 'Phugoid':
    #         tf = t_long * n_periods
    #     elif emotion == 'Shortperiod':
    #         tf = t_short * n_periods
    #     else:
    #         raise TypeError('Eigen motion type not specified')
    # else:
    #     tf *= n_periods
    tf *= n_periods
    t_sim = np.arange(ti, ti + tf, dt)
    deltas = np.append(deltas, np.zeros(len(t_sim)-len(deltas)))
    'First simulation simulates an input pulse taking average deflection, estimated width pulse and is stick fixed'
    yout = response_system(t_sim, sys, "forced", inputvec=deltas)
    vel, alpha, theta, q, delta_e = get_validation(t_sim)
    'Second simulation takes real elevator deflection, thus it is stick free'
    yout2 = response_system(t_sim, sys, "forced", inputvec=delta_e)

    'Sum initial condition'  # integration constant
    initial_condition = np.array([vel[0], alpha[0], theta[0], q[0]]).reshape((4, 1))
    yout += initial_condition
    yout2 += initial_condition
    yout.reshape((4, len(t_sim)))
    yout2.reshape((4, len(t_sim)))
    # print(ysim.shape, t_sim.shape)
    'Readjust plot to show only duration'
    t_plot = tf if t_plot is None else t_plot
    rows = int(t_plot*n_periods/dt)
    t_out = t_sim
    ysim = yout[:, :rows]
    ysim2 = yout2[:, :rows]
    t_sim = t_sim[:rows]
    vel, alpha, theta, q, delta_e = get_validation(t_sim)
    if plot:
        """ Plot all state variables over time """
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(nrows=2, ncols=2, figsize=(9, 7), dpi=100)

        # U vel
        ax1.set_ylabel("Time [s]", fontsize=10)
        ax1.set_xlabel("Velocity [m/s]", fontsize=10)
        ax1.plot(t_sim, ysim[0], label='Simulation')
        ax1.plot(t_sim, ysim2[0], label='Simulation with real $d_{e}$')
        ax1.plot(t_sim, vel, label='Validation')
        ax1.grid(which='major', linestyle=':', color='gray')
        ax1.legend()

        # Alpha
        ax2.set_ylabel("Time [s]", fontsize=10)
        ax2.set_xlabel("Alpha [rad]", fontsize=10)
        ax2.plot(t_sim, ysim[1], label='Simulation')
        ax2.plot(t_sim, ysim2[1], label='Simulation with real $d_{e}$')
        ax2.plot(t_sim, alpha, label='Validation')
        ax2.grid(which='major', linestyle=':', color='gray')
        ax2.yaxis.tick_right()
        ax2.yaxis.set_label_position("right")
        ax2.legend()

        # Pitch angle
        ax3.set_ylabel("Time [s]", fontsize=10)
        ax3.set_xlabel("Theta [rad]", fontsize=10)
        ax3.plot(t_sim, ysim[2], label='Simulation')
        ax3.plot(t_sim, ysim2[2], label='Simulation with real $d_{e}$')
        ax3.plot(t_sim, theta, label='Validation')
        ax3.grid(which='major', linestyle=':', color='gray')
        ax3.legend()

        # Pitch rate
        ax4.set_ylabel("Input deflections [rad]", fontsize=10)
        ax4.set_xlabel("Time elapsed [s]", fontsize=10)
        ax4.plot(t_sim, ysim[3], label='Simulation')
        ax4.plot(t_sim, ysim2[3], label='Simulation with real $d_{e}$')
        ax4.plot(t_sim, q, label='Validation')
        box = ax4.get_position()
        ax4.set_position([box.x0, box.y0 + box.height * 0.22,
                          box.width, box.height * 0.78])
        # Put a legend below current axis
        ax4.legend(loc='upper center', bbox_to_anchor=(0.5, -0.22),
                   fancybox=True, shadow=False, ncol=2)
        ax4.grid(which='major', linestyle=':', color='gray')
        ax4.yaxis.tick_right()
        ax4.yaxis.set_label_position("right")
        plt.suptitle(emotion, fontsize=20)

        plt.savefig(f"plots/{'_'.join(emotion.lower().split(' '))}.png", dpi=600)
        plt.show()
    return t_out, yout, deltas


'Find initial time and duration of eigenmotions'


def find_input_interval(t0, margin, t_sim, delta_sim, emotion):
    t_sim = np.arange(t_sim[0] - margin[0], t_sim[-1]+margin[1], 0.1)
    delta_e = get_validation(t_sim)[4]
    # delta_e -= delta_e[0]
    delta_sim = np.append(np.zeros(int(margin[0]/0.1)-1), delta_sim)
    delta_sim = np.append(delta_sim, np.zeros(int(margin[1]/0.1)))
    plt.plot(t_sim, delta_sim, label="Simulated deflection of the elevator")
    plt.plot(t_sim, delta_e, label="Real deflection of the elevator")
    plt.axvline(x=t0, color="r", label="Estimated start time of " + emotion.lower())
    plt.title("Elevator deflection of " + emotion)
    plt.xlabel("Time [s]", fontsize=10)
    plt.ylabel("$d_e$ [rad]", fontsize=10)
    plt.legend()
    plt.grid()
    plt.savefig(f"plots/elevator_input_{'_'.join(emotion.lower().split(' '))}.png", dpi=600)
    plt.show()
    return

'Initial parameters for symmetric eigenmotion simulation'
'Fetch manually from graphs of the input deflection of elevator the follwing parameters'
delta_ph = -0.006           # rad
tph0 = 3229                 # sec
tf_ph = 22                  # width interval duration of input

delta_shp = -0.06           # rad
tshp0 = 3634                # sec
tf_shp = 0.5 #3645.5-tshp0  # width interval duration of shortperiod/input

'Generate responses'

t_ph, phugoid, sim_ph_delta = sim_sym_resp('Phugoid', tph0, tf_ph, plot=True, n_periods=10, d_e=delta_ph, delay=2)
t_shp, short_period, sim_shp_delta = sim_sym_resp('Shortperiod', tshp0, tf_shp, plot=True, n_periods=30, d_e=delta_shp, delay=0.0)

'Plot both inputs'
# Reduce simulation time of Phugoid as now we are interested only on time during
# input not on the response over time (damping)
t_ph = t_ph[:int(len(t_ph)/4.0)]
sim_ph_delta = sim_ph_delta[:int(len(sim_ph_delta)/4.0)]

find_input_interval(53 * 60 + 57, [25, 0], t_ph, sim_ph_delta, "Phugoid")
find_input_interval(1 * 3600 + 35, [20, 20], t_shp, sim_shp_delta, "Short Period")