import numpy as np


# Auxiliary function to create a meshgrid similar to MATLAB's meshgrid


def meshgrid(x, y):  # good
    x = np.asarray(x).reshape(-1, 1)  # column vector
    y = np.asarray(y).reshape(1, -1)  # row vector
    X = np.repeat(x, y.shape[1], axis=1)
    Y = np.repeat(y, x.shape[0], axis=0)
    return X, Y


def compute_kx(nx, Lx):  # good
    k = np.concatenate((
        np.arange(0, nx // 2 + 1),
        np.arange(-nx // 2 + 1, 0)
    ))
    kx = 1j * k * (2 * np.pi / Lx)
    return kx


def ext(x):  # good

    [nx, ny] = x.shape
    x_ext = np.zeros((2 * nx, 2 * ny))

    # Original block
    x_ext[0:nx, 0:ny] = x

    # Flip horizontally
    x_ext[0:nx, ny:2*ny] = x[:, ::-1]

    # Flip vertically
    x_ext[nx:2*nx, 0:ny] = x[::-1, :]

    # Flip both
    x_ext[nx:2*nx, ny:2*ny] = x[::-1, ::-1]
    return x_ext


def extback(x_ext):  # good
    # Shrinks from 2*nx x 2*ny back to nx x ny (upper-left block)
    [nx_ext, ny_ext] = (x_ext).shape
    nx = int(nx_ext / 2)
    ny = int(ny_ext / 2)
    x_back = x_ext[0:nx, 0:ny]

    return x_back


def f(phi):  # good
    # f = @(x) 0.25*(x.^2-1).^2
    fphi = (phi ** 2 - 1) ** 2 / 4
    return fphi


def f_SAV(phi, gamma0):
    fphi = (phi ** 2 - 1 - gamma0) ** 2 / 4
    return fphi


def r0_fun(phi0, hx, hy, C0, gamma0):

    e1 = np.fft.fft2(f_SAV(phi0, gamma0))
    r0 = np.sqrt(hx * hy * e1[0, 0] + C0)
    # fphi = f(phi0)
    # r0 = np.sqrt(hx * hy * np.sum(fphi) + C0)
    return r0


def df(phi, gamma0):  # good
    return phi ** 3 - (1 + gamma0) * phi


def Lap_SAV(phi, k2, boundary):
    # phi_hat = fft_filtered(phi)
    # result = np.real(np.fft.ifft2(k2 * phi_hat))
    if boundary == "periodic":
        Lphi = np.real(np.fft.ifft2(k2 * np.fft.fft2(phi)))
    elif boundary == "neumann":
        Lphi = np.real(np.fft.ifft2(k2 * fft_filtered(phi)))
    else:
        print("Boundary condition not recognized. Use 'periodic' or 'neumann'.")
    return Lphi

    # return result


def A_inv_CN(phi, dt, k2, k4, gamma0, epsilon2, boundary):
    if boundary == "periodic":
        denom = 1 + ((dt / 2) * epsilon2 * k4 - (dt / 2) * gamma0 * k2)
        Ai = np.real(np.fft.ifft2(np.fft.fft2(phi) / denom))
    elif boundary == "neumann":
        denom = 1 + ((dt / 2) * epsilon2 * k4 - (dt / 2) * gamma0 * k2)
        Ai = np.real(np.fft.ifft2(fft_filtered(phi) / denom))
    return Ai


def fft_filtered(x):  # fft2 is equivalent to fft in julia for 2d array, good
    return np.fft.fft2(x)


def b_fun(phi, hx, hy, C0, gamma0):  # good
    e1 = np.fft.fft2(f_SAV(phi, gamma0))
    return df(phi, gamma0) / np.sqrt(e1[0, 0] * hx * hy + C0)


def g_fun_CN(phi0, r0, b, dt, hx, hy, epsilon2, gamma0, C0, k2, boundary):  # good
    Lap_phi0 = Lap_SAV(phi0, k2, boundary)
    Lap_Lap_phi0 = Lap_SAV(Lap_phi0, k2, boundary)

    bphi0 = np.fft.fft2(b * phi0)
    bphi0 = hx * hy * bphi0[0, 0]

    e1 = np.fft.fft2(f_SAV(phi0, gamma0))

    g = phi0 - (dt / 2) * epsilon2 * Lap_Lap_phi0 + (dt / 2) * gamma0 * Lap_phi0 + dt * Lap_SAV(b, k2, boundary) * \
        (r0 - (1 / 4) * bphi0)

    return g


def r_fun(phi, phi0, r0, b, hx, hy, C0, dt, gamma0):  # good
    bphi0 = np.fft.fft2(b * phi0)
    bphi0 = hx * hy * bphi0[0, 0]

    bphi = np.fft.fft2(b * phi)
    bphi = hx * hy * bphi[0, 0]

    e1 = np.fft.fft2(f_SAV(phi0, gamma0))

    r = r0 + (1 / 2) * (bphi - bphi0)

    return r


def calculate_discrete_energy_sav(phi, h2, epsilon2, k2, gamma0, r, C0):

    E_gamma = h2 * np.fft.fft2(gamma0 / 2 * phi ** 2)
    E_gamma = E_gamma[1, 1]

    a = r ** 2 - C0 - (gamma0 ** 2 + 2 * gamma0) / 4

    a = h2 * np.sum(f_SAV(phi, gamma0))
    sum_i = np.sum((phi[1:, :] - phi[0:-1, :]) ** 2)

    b = (epsilon2 / 2) * sum_i
    sum_j = np.sum((phi[:, 2:] - phi[:, 1:-1]) ** 2)

    c = (epsilon2 / 2) * sum_j
    E = a + b + c + E_gamma

    return E


def ch_r_error(r, phi, h2, C0, gamma0):
    e1 = np.fft.fft2(f_SAV(phi, gamma0))
    D = r - np.sqrt(h2 * e1[0, 0] + C0)
    D = D / np.sqrt(h2 * e1[0, 0] + C0)
    return D
