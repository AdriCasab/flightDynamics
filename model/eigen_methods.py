"""
Group 1203F1
Flight Dynamics Assignment - Useful methods for eigenmotions
"""

from math import *
import numpy as np
import pandas as pd


def compute_eigenmotion_parameters(eigenvalue):
    """
    Compute the eigenmotion stability parameters associated to an eigenvalue.
    Input should be a complex eigenvalue
    """
    omega_n = sqrt(eigenvalue.real ** 2 + eigenvalue.imag ** 2)
    damp_ratio = -eigenvalue.real / sqrt(eigenvalue.real ** 2 + eigenvalue.imag ** 2)

    if eigenvalue.real > 0:
        # eigenmotion is divergent
        t_double = log(2) / eigenvalue.real
        convergence = 'Divergent'

        if np.isclose(eigenvalue.imag, 0.):
            time_constant = 1. / ((damp_ratio - sqrt(damp_ratio ** 2 - 1)) * omega_n)
            period = None
            periodicity = 'Aperiodic'
        else:
            period = float(2*pi / eigenvalue.imag)
            time_constant = None
            periodicity = 'Periodic'

        return period, damp_ratio, t_double, omega_n, time_constant, convergence, periodicity

    elif eigenvalue.real < 0:
        t_half = log(0.5) / eigenvalue.real
        convergence = 'Convergent'

        if np.isclose(eigenvalue.imag, 0.):
            time_constant = 1. / ((damp_ratio - sqrt(damp_ratio ** 2 - 1)) * omega_n)
            period = None
            periodicity = 'Aperiodic'
        else:
            period = float(2 * pi / eigenvalue.imag)
            time_constant = None
            periodicity = 'Periodic'

        return period, damp_ratio, t_half, omega_n, time_constant, convergence, periodicity


def characterize_eigenmotions(A):
    """
    :param A: 4x4 state matrix
    :return: Dataframe with eigenmotions and corresponding characteristics
    """
    eigenvals = np.linalg.eig(A)[0]
    eigen_pars = np.array([compute_eigenmotion_parameters(x) for x in eigenvals])

    # Generate dataframe with eigenmode characteristics
    eigen_df = pd.DataFrame()
    eigen_df['Eigenmotion_nr'] = list(range(1, 5))
    eigen_df['Real'] = eigenvals.real
    eigen_df['Imag'] = eigenvals.imag
    eigen_df['Convergence'] = eigen_pars[:, 5]
    eigen_df['Periodicity'] = eigen_pars[:, 6]
    eigen_df['Period'] = eigen_pars[:, 0]
    eigen_df['Damp_ratio'] = eigen_pars[:, 1]
    eigen_df['T_half/double'] = eigen_pars[:, 2]
    eigen_df['Omega_n'] = eigen_pars[:, 3]
    eigen_df['Time_constant'] = eigen_pars[:, 4]
    eigen_df.set_index('Eigenmotion_nr', drop=True, inplace=True)

    return eigen_df


