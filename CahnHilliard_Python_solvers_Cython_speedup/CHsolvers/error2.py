from . import aux_functions_NMG as aux
import numpy as np

from .laplacian import laplace


def error2(c_old, c_new, mu, nxt, nyt, dt, xright, xleft, yright, yleft, boundary):
    """
    Calculate the residual for phi
    :param c_old: old phi
    :param c_new: updated phi
    :param mu: updated mu
    :param nxt: temp (number of grid points in x-direction, locally defined)
    :param nyt: temp (number of grid points in y-direction, locally defined)
    :return: res2, Frobenius norm (residual), calculated after each vcycle update to c_new and mu
    """
    rr = aux.dmatrix(nxt, nyt)
    for i in range(nxt):
        for j in range(nyt):
            rr[i][j] = mu[i][j] - c_old[i][j]

    sor = laplace(rr, nxt, nyt, xright, xleft, yright, yleft, boundary)
    for i in range(nxt):
        for j in range(nyt):
            rr[i][j] = sor[i][j] - (c_new[i][j] - c_old[i][j]) / dt
    res2 = np.sqrt(np.sum(rr**2) / (nxt * nyt))
    return res2
