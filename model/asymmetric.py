"""
Group 1203F1 Flight Dynamics Assignment
Asymmetric Numerical Model - Characterizing and Simulating Eigenmotions
"""

from numpy.linalg import inv, eig
from Cit_par import *
import matplotlib.pyplot as plt
import numpy as np
from control import ss, forced_response
from datetime import datetime
import os
from eigen_methods import characterize_eigenmotions

# Directory in which response plots are to be saved
directory = "plots/"

if not os.path.exists(directory):
    os.makedirs(directory)


def get_asymmetric_state_space(time):
    '''
    from given time return state space for asymmetric motion
    :param time: import time in seconds
    :return: matrix A and state space for asymmetric motion
    '''
    V, mub, CL, muc, CX0, CZ0 = get_stationary_condition(time)

    P = np.array([[(CYbdot-2.*mub) * b / V, 0, 0, 0],
                  [0, -0.5 * b / V, 0, 0],
                  [0, 0, -2. * mub * KX2 * (b / V) ** 2., 2. * mub * KXZ * (b / V) ** 2.],
                  [Cnbdot * b / V, 0, 2. * mub * KXZ * (b / V) ** 2., -2. * mub * KZ2 * (b / V) ** 2.]])

    Q = np.array([[-CYb, -CL, -CYp * b / (2. * V), -b * (CYr - 4 * mub) / (2. * V)],
                  [0, 0, -b / (2. * V), 0],
                  [-Clb, 0, -Clp * b / (2. * V), -Clr * b / (2 * V)],
                  [-Cnb, 0, -Cnp * b / (2. * V), -Cnr * b / (2. * V)]])

    R = np.array([[-CYda, -CYdr],
                   [0, 0],
                   [-Clda, -Cldr],
                   [-Cnda, -Cndr]])

    # Set up state-space system

    A = inv(P) @ Q                  # state matrix
    B = inv(P) @ R                  # input matrix
    C = np.identity(A.shape[0])     # output matrix
    D = np.zeros(B.shape)           # direct matrix

    sys = ss(A, B, C, D)            # assemble state-space system using control module

    return A, sys


def make_pulse_input(magnitude, start_time, end_time, t_vector):
    """ Generates square pulse signal active between start_time and end_time
    """
    return [magnitude*(start_time < entry < end_time) for entry in t_vector]


def make_step_input(magnitude, start_time, t_vector):
    """ Generates step signal starting from start_time
    """
    return [magnitude * (start_time < entry) for entry in t_vector]


""" 
Simulation methods for verification purposes -- the user can specify properties of the pulse or step inputs
to be generated.
Useful to verify whether state-space system responds as expected for certain inputs and eigenmodes.
"""


def verification_sim_dutch_roll(sys, pulse_mag, pulse_start, pulse_end, sim_start, sim_end, step, plot=True):
    """
    Simulate dutch roll eigenmotion:
        --> input vector is impulse rudder deflection of [pulse_mag] active between [pulse_start] and [pulse_end]
        --> response is simulated between [sim_start] and [sim_end] with a stepsize of [step]
        --> NOTE: make sure pulse times match sim times (on same time-scale)
    """
    # Simulation time vector
    t_sim = np.arange(sim_start, sim_end, step)
    # Assemble input vector
    input_vector = np.vstack((np.zeros(len(t_sim)), make_pulse_input(pulse_mag, pulse_start, pulse_end, t_sim)))
    # Simulate response
    t, yout, xout = forced_response(sys, t_sim, U=input_vector)

    if plot:
        """ Plot all state variables over time """
        # Sideslip angle (beta)
        plt.subplot(221)
        plt.plot(t_sim, yout[0])
        plt.title("Beta vs t")
        # Roll angle (phi)
        plt.subplot(222)
        plt.plot(t_sim, yout[1])
        plt.title("Phi vs t")
        # Roll rate (p)
        plt.subplot(223)
        plt.plot(t_sim, yout[2])
        plt.title("Roll rate vs t")
        # Yaw rate (r)
        plt.subplot(224)
        plt.plot(t_sim, yout[3])
        plt.title("Yaw rate vs t")

    return yout


def verification_sim_aperiodic_roll(sys, step_mag, step_start, sim_start, sim_end, step, plot=True):
    """
        Simulate aperiodic roll eigenmotion:
            --> input vector is step aileron deflection of [step_mag] active from [step_start]
            --> response is simulated between [sim_start] and [sim_end] with a stepsize of [step]
            --> NOTE: make sure step times match sim times (on same time-scale)
    """
    # Simulation time vector
    t_sim = np.arange(sim_start, sim_end, step)
    # Assemble input vector
    input_vector = np.vstack((make_step_input(step_mag, step_start, t_sim), np.zeros(len(t_sim))))
    # Simulate response
    t, yout, xout = forced_response(sys, t_sim, U=input_vector)

    if plot:
        """ Plot all state variables over time """
        # Sideslip angle (beta)
        plt.subplot(221)
        plt.plot(t_sim, yout[0])
        plt.title("Beta vs t")
        # Roll angle (phi)
        plt.subplot(222)
        plt.plot(t_sim, yout[1])
        plt.title("Phi vs t")
        # Roll rate (p)
        plt.subplot(223)
        plt.plot(t_sim, yout[2])
        plt.title("Roll rate vs t")
        # Yaw rate (r)
        plt.subplot(224)
        plt.plot(t_sim, yout[3])
        plt.title("Yaw rate vs t")


def verification_sim_spiral(sys, pulse_mag, pulse_start, pulse_end, sim_start, sim_end, step, plot=True):
    """
        Simulate spiral eigenmotion (unstable):
            --> input vector is impulse aileron deflection of [pulse_mag] active between [pulse_start] and [pulse_end]
            --> response is simulated between [sim_start] and [sim_end] with a stepsize of [step]
            --> NOTE: make sure pulse times match sim times (on same time-scale)
        """
    # Simulation time vector
    t_sim = np.arange(sim_start, sim_end, step)
    # Assemble input vector
    input_vector = np.vstack((make_pulse_input(pulse_mag, pulse_start, pulse_end, t_sim), np.zeros(len(t_sim))))
    # Simulate response
    t, yout, xout = forced_response(sys, t_sim, U=input_vector)

    if plot:
        """ Plot all state variables over time """
        # Sideslip angle (beta)
        plt.subplot(221)
        plt.plot(t_sim, yout[0])
        plt.title("Beta vs t")
        # Roll angle (phi)
        plt.subplot(222)
        plt.plot(t_sim, yout[1])
        plt.title("Phi vs t")
        # Roll rate (p)
        plt.subplot(223)
        plt.plot(t_sim, yout[2])
        plt.title("Roll rate vs t")
        # Yaw rate (r)
        plt.subplot(224)
        plt.plot(t_sim, yout[3])
        plt.title("Yaw rate vs t")

    return yout


""" 
Simulation methods for numerical model -- these use input vectors taken from the Flight Test Data. 
"""


def sim_arbitrary_input(sys, input_vector, t_sim, ic=None):
    """ Returns [beta, phi, roll rate, yaw rate]"""
    if ic is not None:
        _, yout  = forced_response(sys, t_sim, U=input_vector, X0=ic)
    else:
        _, yout = forced_response(sys, t_sim, U=input_vector)
    return yout


def sim_asymmetric_eigenmotion(start_time, end_time, motion='Dutch Roll', step=0.1, plot=True):
    """ Simulate  eigenmotion using inputs from Flight Test Data"""

    t_sim = np.arange(start_time, end_time, step)
    drs, das, roll_angles, roll_rates, yaw_rates = [], [], [], [], []

    for time in t_sim:
        drs.append(get_arbitrary_variable(time, d_r))
        das.append(get_arbitrary_variable(time, d_a))
        roll_angles.append(get_arbitrary_variable(time, roll_angle))
        roll_rates.append(get_arbitrary_variable(time, roll_rate))
        yaw_rates.append(get_arbitrary_variable(time, yaw_rate))

    A, sys = get_asymmetric_state_space(t_sim[0])
    print(characterize_eigenmotions(A))

    if motion == 'Dutch Roll':
        """ 
        Use rudder deflection impulse from flight test data - assume no aileron deflection for simulation (this will
        be source of deviation from validation, as in reality aileron deflection will not be perfectly zero.
        """
        # calibration_factor = 1.2215
        # u = np.vstack((np.array(das) - das[0], -(np.array(drs) - calibration_factor * drs[0])))
        # u = np.vstack((np.array(das) - das[0], -(np.array(drs))))
        u = np.vstack((np.array(das) - das[0], -(np.array(drs)-drs[0])))
        yout = sim_arbitrary_input(sys, u, t_sim)

    elif motion == 'Aperiodic Roll':
        """
        Use aileron deflection step from flight test data- assume no rudder deflection (see previous docstring for note)
        """
        # u = np.vstack((np.array(das)-das[0], np.zeros(len(t_sim))))
        u = np.vstack((np.array(das)-das[0], np.zeros(len(t_sim))))
        yout = sim_arbitrary_input(sys, u, t_sim, ic=[0., roll_angles[0], roll_rates[0], yaw_rates[0]])

    elif motion == 'Spiral':
        """
        Use aileron deflection impulse from flight test data - assume no rudder deflection 
        (see previous docstring for note)
        """
        # u = np.vstack((np.array(das)-das[0], np.zeros(len(t_sim))))
        # u = np.vstack((np.array(das)-das[0]/1.1, np.array(drs)-drs[0]))

        calibration_factor = 1.2215
        u = np.vstack((np.array(das) - das[0], -(np.array(drs)-calibration_factor*drs[0])))
        yout = sim_arbitrary_input(sys, u, t_sim, ic=[0., roll_angles[0],  roll_rates[0], yaw_rates[0]])

    else:
        u = None
        raise ValueError("Invalid eigenmotion -- par [motion] must be 'Dutch Roll', 'Aperiodic Roll', or 'Spiral'")

    if plot:
        # Plot state variables
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(nrows=2, ncols=2, figsize=(9, 7), dpi=100)
        ax1.set_ylabel("Roll angle [rad]", fontsize=10)
        ax1.set_xlabel("Time elapsed [s]", fontsize=10)
        ax1.plot(np.array(t_sim)-t_sim[0], roll_angles, label='Validation')
        ax1.plot(np.array(t_sim)-t_sim[0], yout[1], label='Simulation')
        ax1.grid(which='major', linestyle=':', color='gray')
        ax1.set_xlim(xmin=0)
        ax1.legend()

        ax2.set_ylabel("Roll rate [rad/s]", fontsize=10)
        ax2.set_xlabel("Time elapsed [s]", fontsize=10)
        ax2.plot(np.array(t_sim) - t_sim[0], roll_rates, label='Validation')
        ax2.plot(np.array(t_sim) - t_sim[0], yout[2], label='Simulation')
        ax2.grid(which='major', linestyle=':', color='gray')
        ax2.yaxis.tick_right()
        ax2.yaxis.set_label_position("right")
        ax2.set_xlim(xmin=0)
        ax2.legend()

        ax3.set_ylabel("Yaw rate [rad/s]", fontsize=10)
        ax3.set_xlabel("Time elapsed [s]", fontsize=10)
        ax3.plot(np.array(t_sim) - t_sim[0], yaw_rates, label='Validation')
        ax3.plot(np.array(t_sim) - t_sim[0], yout[3], label='Simulation')
        ax3.grid(which='major', linestyle=':', color='gray')
        ax3.set_xlim(xmin=0)
        ax3.legend()

        ax4.set_ylabel("Input deflections [rad]", fontsize=10)
        ax4.set_xlabel("Time elapsed [s]", fontsize=10)
        ax4.plot(np.array(t_sim) - t_sim[0], drs, label='Measured rudder')
        if motion == 'Spiral':
            ax4.plot(np.array(t_sim) - t_sim[0], np.array(drs) - calibration_factor*drs[0], label='Calibrated rudder')
        else:
            ax4.plot(np.array(t_sim) - t_sim[0], np.array(drs) - drs[0], label='Calibrated rudder')
        ax4.plot(np.array(t_sim) - t_sim[0], das, label='Measured aileron')
        ax4.plot(np.array(t_sim) - t_sim[0], np.array(das) - das[0], label='Calibrated aileron')
        box = ax4.get_position()
        ax4.set_position([box.x0, box.y0 + box.height * 0.22,
                         box.width, box.height * 0.78])
        # Put a legend below current axis
        ax4.legend(loc='upper center', bbox_to_anchor=(0.5, -0.22),
                   fancybox=True, shadow=False, ncol=2)
        ax4.grid(which='major', linestyle=':', color='gray')
        ax4.yaxis.tick_right()
        ax4.yaxis.set_label_position("right")
        ax4.set_xlim(xmin=0)
        plt.suptitle(motion, fontsize=20)

        plt.savefig(f"plots/{'_'.join(motion.lower().split(' '))}.png", dpi=600)
        plt.show()

    return yout


def get_seconds(timestring):
    """ timestring should be in form 'HH:MM:SS' eg: """
    secs = [3600, 60, 1]
    return sum([a * b for a, b in zip(secs, map(int, timestring.split(':')))])


# print(characterize_eigenmotions(get_asymmetric_state_space(3000)[0]))

sim_asymmetric_eigenmotion(get_seconds('01:01:57'), get_seconds('01:02:10'), motion='Dutch Roll', plot=True)
sim_asymmetric_eigenmotion(get_seconds('00:59:10'), get_seconds('00:59:23'), motion='Aperiodic Roll', plot=True)
sim_asymmetric_eigenmotion(get_seconds('01:05:00'), get_seconds('01:07:00'), motion='Spiral', plot=True)
