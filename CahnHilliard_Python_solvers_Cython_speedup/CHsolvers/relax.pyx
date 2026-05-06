import numpy as np
cimport numpy as np


def relax(double[:, :] c_new, double[:, :] mu_new, double[:, :] su, double[:, :] sw,
          int nxt, int nyt, int c_relax,
          double xright, double xleft, double yright, double yleft,
          double dt, double epsilon2, str boundary):
    cdef int iter, i, j
    cdef int periodic = (boundary == "periodic")
    cdef double ht2, x_fac, y_fac, det
    cdef double a0, a1, a2, a3
    cdef double f0, f1

    ht2 = ((xright - xleft) / nxt) ** 2   # ← now comes AFTER all cdefs

    for iter in range(c_relax):

        for i in range(nxt):
            for j in range(nyt):
                if boundary == "neumann":
                    if i > 0 and i < nxt - 1:
                        x_fac = 2.0
                    else:
                        x_fac = 1.0
                    if j > 0 and j < nyt - 1:
                        y_fac = 2.0
                    else:
                        y_fac = 1.0
                elif boundary == "periodic":
                    x_fac = 2.0
                    y_fac = 2.0

                a0 = 1 / dt
                a1 = (x_fac + y_fac) / ht2
                a2 = -(x_fac + y_fac) * epsilon2 / \
                    ht2 - 3 * (c_new[i][j]) ** 2
                a3 = 1.0

                f0 = su[i, j]
                f1 = sw[i, j] - 2 * (
                    c_new[i, j] ** 3
                )  # replaced from c code with a more condensed version
                # boundary cases are slightly different because i-1 doesn't exist for i = 0, for example (same for j)
                if i > 0:
                    f0 += mu_new[i - 1, j] / ht2
                    f1 -= epsilon2 * c_new[i - 1, j] / ht2
                elif boundary == "periodic":
                    # should stay nxt-1 because python is 0-indexed
                    f0 += mu_new[nxt-1, j] / ht2
                    f1 -= epsilon2 * c_new[nxt - 1, j] / ht2
                if i < nxt - 1:
                    f0 += mu_new[i + 1, j] / ht2
                    f1 -= epsilon2 * c_new[i + 1, j] / ht2
                elif boundary == "periodic":
                    f0 += mu_new[0, j] / ht2  # changed
                    f1 -= epsilon2 * c_new[0, j] / ht2  # changed
                if j > 0:
                    f0 += mu_new[i, j - 1] / ht2
                    f1 -= epsilon2 * c_new[i, j - 1] / ht2
                elif boundary == "periodic":
                    f0 += mu_new[i, nyt - 1] / ht2
                    f1 -= epsilon2 * c_new[i, nyt - 1] / ht2
                if j < nyt - 1:
                    f0 += mu_new[i, j + 1] / ht2
                    f1 -= epsilon2 * c_new[i, j + 1] / ht2
                elif boundary == "periodic":
                    f0 += mu_new[i, 0] / ht2  # changed
                    f1 -= epsilon2 * c_new[i, 0] / ht2  # changed
                det = a0 * a3 - a1 * a2
                c_new[i, j] = (a3 * f0 - a1 * f1) / det
                mu_new[i, j] = (-a2 * f0 + a0 * f1) / det
    return c_new, mu_new
