import numpy as np
from . import NMG_solver as solver
import os
from .aux_functions_NMG import *


@time_and_mem
def CahnHilliard_NMG(phi0, t_iter=1e3, dt=2.5e-5, solver_iter=1e4, tol=1e-5, dt_out=10, m=8, epsilon2=np.nan, boundary="periodic", c_relax=2, domain=[1, 0, 1, 0], printres=False, printphi=False, pathname="cd"):
    nx, ny = phi0.shape
    xright, xleft, yright, yleft = domain

    n_level = int(np.log(nx) / np.log(2.0) + 0.1)
    h2 = ((xright - xleft) / nx) * ((yright - yleft) / ny)  # Define mesh size
    if np.isnan(epsilon2):
        # Define ϵ^2 if not prespecified
        epsilon2 = h2 * m ** 2 / (2 * np.sqrt(2) * np.arctanh(0.9)) ** 2
    else:
        # Else overwrite m
        m = np.sqrt((epsilon2 * (2 * np.sqrt(2) * np.arctanh(0.9)) ** 2) / h2)

    mu = np.zeros((nx, ny))  # µ defined as a global variable
    phi_new = phi0.copy()
    phi_old = phi0.copy()

    # Logical index for the need to downsample
    downsampled = nx * ny * t_iter / dt_out > 1e9
    n_timesteps = np.floor(t_iter / dt_out)

    # Save or print initial phi
    if printphi:
        mass_t = np.zeros(int(n_timesteps+1))
        E_t = np.zeros(int(n_timesteps+1))
        t_out = np.arange(0, (int(n_timesteps+1))*dt*dt_out, dt_out*dt)
        if pathname == "cd":
            pathname = os.pwd()
        with open(f"{pathname}phi.csv", "w") as f:
            for i in range(nx):
                for j in range(ny):
                    f.write(f"{phi0[i][j]},")
                f.write("\n")  # end of each row
        phi_t = np.nan
    else:
        if downsampled:
            # we need to round up to ensure we have enough space
            new_dt_out = np.ceil(nx * ny * t_iter / 1e9)
            print("Variable phi_t is too large with dt_out = %4.0f. Downsampling to every %4.0f time steps\n", dt_out, new_dt_out)
            dt_out = new_dt_out
            n_timesteps = np.floor(t_iter / dt_out)
        mass_t = np.zeros(int(n_timesteps+1))
        E_t = np.zeros(int(n_timesteps+1))
        t_out = np.arange(0, (int(n_timesteps+1))*dt*dt_out, dt_out*dt)
        phi_t = np.zeros((nx, ny, int(n_timesteps+1)))
        phi_t[:, :, 0] = phi0

    mass_t[0] = calculate_mass(phi0, h2, nx, ny)
    E_t[0] = calculate_discrete_energy(phi0, h2, epsilon2)

    if printres:
        print("Saving squared residuals per iteration to file in the output directory\n")
    rr = np.zeros((nx, ny))
    for it in range(t_iter):
        print(it)
        solver.cahn(phi_old, phi_new, mu, nx, ny, dt, solver_iter, tol, c_relax, xright,
                    xleft, yright, yleft, epsilon2, n_level, boundary, printres=printres, pathname=pathname)
        phi_old = phi_new.copy()
        if it / t_iter * 100 % 5 == 0:
            print((f"{it / t_iter * 100:3.0f} percent complete"))
        if it % dt_out == 0:
            t_index = int(np.floor(it / dt_out) + 1)
            if printphi:
                with open(f"{pathname}phi.csv", "a") as f:
                    for i in range(nx):
                        for j in range(ny):
                            f.write(f"{phi_new[i][j]},")
                        f.write("\n")  # end of each row
            else:
                phi_t[:, :, t_index] = phi_new
            mass_t[t_index] = calculate_mass(phi_new, h2, nx, ny)
            E_t[t_index] = calculate_discrete_energy(
                phi_new, h2,  epsilon2)
    delta_mass_t = mass_t - mass_t[0]
    E_t = E_t / E_t[0]
    return t_out, phi_t, delta_mass_t, E_t
