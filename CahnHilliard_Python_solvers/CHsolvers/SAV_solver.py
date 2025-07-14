import numpy as np
from . import aux_functions_SAV as aux   # Import auxiliary functions

# This function uses the sav method to solve the
# Cahn-Hilliard equation for the next time step.

# INPUTS
# phi_old  = Prior chemical state.
# r_old    = Prior SAV parameter.
# hx, hy   = Grid spacings in x and y directions.
# k2, k4   = Laplacian and biharmonic operators.
# dt       = Time step.
# epsilon2 = Square of the interface width parameter.
# boundary = 'periodic' or 'neumann' (not explicitly used in this snippet).
# C0       = Regularization parameter.
# gamma0   = Stabilization parameter.

# OUTPUT
# phi_new = Next chemical state.
# r_new   = Next sav state.


def sav_solver(phi_old, phi_prev, r_old, hx, hy, k2, k4, dt, epsilon2, boundary, C0, gamma0, eta, xi_flag, i):

    phi0 = phi_old
    r0 = r_old

    phi0_df = aux.df(phi0, gamma0)  # df at phi0
    Lap_dfphi0 = aux.Lap_SAV(phi0_df, k2, boundary)  # Lap of df(phi0)

    if i == 0:
        phi_bar = aux.A_inv_CN(phi0 + dt/2 * Lap_dfphi0,
                               dt, k2, k4, gamma0, epsilon2, boundary)
    elif i >= 1:
        phi_bar = 1.5 * phi_old - 0.5 * phi_prev
        phi_bar = np.maximum(-1, np.minimum(1, phi_bar))

    # Step 1
    b = aux.b_fun(phi_bar, hx, hy, C0, gamma0)
    g = aux.g_fun_CN(phi0, r0, b, dt, hx, hy, epsilon2,
                     gamma0, C0, k2, boundary)

    AiLb = aux.A_inv_CN(aux.Lap_SAV(b, k2, boundary),
                        dt, k2, k4, gamma0, epsilon2, boundary)
    Aig = aux.A_inv_CN(g, dt, k2, k4, gamma0, epsilon2, boundary)

    gamma = -np.fft.fft2(b*AiLb)
    gamma = gamma[0, 0]*hx*hy

    # Step 2
    bphi = np.fft.fft2(b*Aig)
    bphi = bphi[0, 0]*hx*hy/(1+dt/4*gamma)

    # Step 3
    phi_new = dt/4*bphi*AiLb + Aig
    r_new = aux.r_fun(phi_new, phi_old, r0, b, hx, hy, C0, dt, gamma0)

    # calculate a,b,c
    # Q_phi_new
    E1_new = np.fft.fft2(aux.f_SAV(phi_new, gamma0))
    E1_new = E1_new[0, 0]*hx*hy
    Q_phi_new = np.sqrt(E1_new+C0)
    # q_tilde
    q_tilde = r_new
    # mu_half and Lap_mu_half
    phi_half = (phi_new + phi_old) / 2
    r_half = (r_new + r_old) / 2
    mu_half = aux.Lap_SAV(phi_half, k2, boundary) + r_half * b
    Lap_mu_half = aux.Lap_SAV(mu_half, k2, boundary)
    # muGmu_half
    muGmu_half = np.fft.fft2(mu_half * (-Lap_mu_half))
    muGmu_half = muGmu_half[0, 0]*hx*hy
    # a, b, c
    a = (q_tilde - Q_phi_new) ** 2
    b = 2 * (q_tilde - Q_phi_new) * Q_phi_new
    c = Q_phi_new ** 2 - q_tilde ** 2 - dt * eta * muGmu_half

    if xi_flag == 0:
        xi = 1  # No relaxation when xi_flag is 1
    elif xi_flag == 1:
        if a > 0:
            discriminant = b**2 - 4 * a * c  # Calculate discriminant
            if discriminant >= 0:
                xi = np.real((-b - np.sqrt(discriminant)) / (2 * a))
                # Restrict xi to the range [0,1]
                xi = np.max([0, np.min([xi, 1])])
            else:
                xi = 1  # When discriminant < 0, choose xi = 1
        else:
            # Special case when a = 0
            if b != 0:
                xi = -c / b
                # Restrict xi to the range [0,1]
                xi = np.max([0, np.min([xi, 1])])
            else:
                if c <= 0:
                    xi = 1  # When c <= 0, choose xi = 1
                else:
                    xi = 0  # When c > 0, no solution, choose xi = 0

    r_new = xi * q_tilde + (1 - xi) * Q_phi_new

    return phi_new, r_new
