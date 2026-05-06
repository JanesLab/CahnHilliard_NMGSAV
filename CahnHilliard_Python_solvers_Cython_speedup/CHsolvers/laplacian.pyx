from . import aux_functions_NMG as aux

import numpy as np
cimport numpy as np


def laplace(double[:, :] a, int nxt, int nyt,
            double xright, double xleft, double yright, double yleft,
            str boundary):
    cdef int i, j
    cdef int periodic = (boundary == "periodic")
    cdef double h2, dadx_L, dadx_R, dady_B, dady_T
    cdef double[:, :] lap_a = np.empty((nxt, nyt), dtype=np.float64)
    h2 = ((xright - xleft) / nxt) ** 2
    for i in range(nxt):
        for j in range(nyt):
            if i > 0:
                dadx_L = a[i, j] - a[i - 1, j]
            else:
                if boundary == "neumann":
                    dadx_L = 0
                elif boundary == "periodic":
                    dadx_L = a[i, j] - a[nxt-1, j]
            if i < nxt - 1:
                dadx_R = a[i + 1, j] - a[i, j]
            else:
                if boundary == "neumann":
                    dadx_R = 0
                elif boundary == "periodic":
                    dadx_R = a[0, j] - a[i, j]
            if j > 0:
                dady_B = a[i, j] - a[i, j - 1]
            else:
                if boundary == "neumann":
                    dady_B = 0
                elif boundary == "periodic":
                    dady_B = a[i, j] - a[i, nyt-1]
            if j < nyt - 1:
                dady_T = a[i, j + 1] - a[i, j]
            else:
                if boundary == "neumann":
                    dady_T = 0
                elif boundary == "periodic":
                    dady_T = a[i, 0] - a[i, j]
            lap_a[i, j] = (dadx_R - dadx_L + dady_T - dady_B) / h2
    return lap_a
