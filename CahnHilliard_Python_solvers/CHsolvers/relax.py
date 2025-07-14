import numpy as np


def relax(
    c_new, mu_new, su, sw, nxt, nyt, c_relax, xright, xleft, yright, yleft, dt, epsilon2, boundary,
):
    """
    SMOOTH Relaxation operator. This is just solving x =b*A-1 for the system of equations c_new and mu_new, where A is
    the LHS of equations 22 and 23, and b is the RHS.
    :param c_new: c to be smoothed
    :param mu_new: mu to be smoothed
    :param su: sc, locally defined
    :param sw: smu, locally defined
    :param nxt: temp (number of grid points in x-direction, locally defined)
    :param nyt: temp (number of grid points in y-direction, locally defined)
    :param c_relax: number of relaxation operations
    :param xright: right x-coordinate
    :param xleft: left x-coordinate
    :param dt: time step
    :param epsilon2: ϵ^2
    :return: c_new, mu_new
    """
    ht2 = ((xright - xleft) / nxt) ** 2  # h2 temp, defined locally
    a = np.empty(4)
    f = np.empty(2)
    # print("c_new before relaxation: \n", c_new)
    # print("mu_new before relaxation: \n", mu_new)
    # print("su before relaxation: \n", su)
    # print("sw before relaxation: \n", sw)
    for iter in range(c_relax):  # c_relax is defined to be 2 in CHsolver.c
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

                a[0] = 1 / dt
                a[1] = (x_fac + y_fac) / ht2
                a[2] = -(x_fac + y_fac) * epsilon2 / \
                    ht2 - 3 * (c_new[i][j]) ** 2
                a[3] = 1.0

                f[0] = su[i][j]
                f[1] = sw[i][j] - 2 * (
                    c_new[i][j] ** 3
                )  # replaced from c code with a more condensed version
                # boundary cases are slightly different because i-1 doesn't exist for i = 0, for example (same for j)
                if i > 0:
                    f[0] += mu_new[i - 1][j] / ht2
                    f[1] -= epsilon2 * c_new[i - 1][j] / ht2
                elif boundary == "periodic":
                    # should stay nxt-1 because python is 0-indexed
                    f[0] += mu_new[nxt-1][j] / ht2
                    f[1] -= epsilon2 * c_new[nxt - 1][j] / ht2
                if i < nxt - 1:
                    f[0] += mu_new[i + 1][j] / ht2
                    f[1] -= epsilon2 * c_new[i + 1][j] / ht2
                elif boundary == "periodic":
                    f[0] += mu_new[0][j] / ht2  # changed
                    f[1] -= epsilon2 * c_new[0][j] / ht2  # changed
                if j > 0:
                    f[0] += mu_new[i][j - 1] / ht2
                    f[1] -= epsilon2 * c_new[i][j - 1] / ht2
                elif boundary == "periodic":
                    f[0] += mu_new[i][nyt - 1] / ht2
                    f[1] -= epsilon2 * c_new[i][nyt - 1] / ht2
                if j < nyt - 1:
                    f[0] += mu_new[i][j + 1] / ht2
                    f[1] -= epsilon2 * c_new[i][j + 1] / ht2
                elif boundary == "periodic":
                    f[0] += mu_new[i][0] / ht2  # changed
                    f[1] -= epsilon2 * c_new[i][0] / ht2  # changed
                det = a[0] * a[3] - a[1] * a[2]
                c_new[i][j] = (a[3] * f[0] - a[1] * f[1]) / det
                mu_new[i][j] = (-a[2] * f[0] + a[0] * f[1]) / det
    return c_new, mu_new
