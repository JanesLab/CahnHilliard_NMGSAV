from . import aux_functions_NMG as aux


def laplace(a, nxt, nyt, xright, xleft, yright, yleft, boundary):
    """
    Compute the discrete Laplacian of a
    :param xright:
    :param a: matrix
    :param nxt: nx temp (number of grid points in x-direction, locally defined)
    :param nyt: ny temp (number of grid points in y-direction, locally defined)
    :return: lap_a, the discrete laplacian of a
    """
    lap_a = aux.dmatrix(nxt, nyt)
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
