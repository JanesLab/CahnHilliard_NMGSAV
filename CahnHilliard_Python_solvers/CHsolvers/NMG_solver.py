import numpy as np
from . import aux_functions_NMG as aux
from . import laplacian as ll
from . import error2
from . import relax

# note: the C code updates ct, but it doesn't appear to be used anywhere. We'll need to check if it should be returned.


def source(c_old, nx, ny, dt, xright, xleft, yright, yleft, boundary):
    """
    Compute the source term for phi and mu
    :param c_old: phi at a time step
    :return: src_c, the source term for phi, and src_mu, the source term for mu
    """

    src_mu = np.zeros((nx, ny))
    src_c = aux.dmatrix(nx, ny)
    ct = ll.laplace(c_old, nx, ny, xright, xleft, yright, yleft, boundary)
    src_c = c_old / dt - ct  # update source term of phi
    # for i in range(nx):
    #     for j in range(ny):
    #         src_c[i, j] = c_old[i, j] / dt - \
    #             ct[i, j]  # update source term of phi
    #         src_mu[i, j] = 0  # set source term for mu to zero
    return src_c, src_mu


def restrict_ch(uf, vf, nxc, nyc):
    """
    Restrict the defect twofold in each direction
    uf and vf get compressed to uc and vc with dimensions nxc and nyc
    Note that changing from C to Python requires adding 1 instead of subtracting in formulas
    :param uf: uf matrix to be restricted
    :param vf: vf matrix to be restricted
    :param nxc: number of grid points in x-direction of uc
    :param nyc: number of grid points in y-direction of vc
    :return: uc, vc
    """
    uc = aux.dmatrix(nxc, nyc)
    vc = aux.dmatrix(nxc, nyc)
    uc = 0.25 * (uf[::2, ::2] + uf[1::2, ::2] + uf[::2, 1::2] + uf[1::2, 1::2])
    vc = 0.25 * (vf[::2, ::2] + vf[1::2, ::2] + vf[::2, 1::2] + vf[1::2, 1::2])
    # for i in range(nxc):
    #     for j in range(nyc):
    #         uc[i][j] = 0.25 * (
    #             uf[2 * i][2 * j]
    #             + uf[2 * i + 1][2 * j]
    #             + uf[2 * i][2 * j + 1]
    #             + uf[2 * i + 1][2 * j + 1]
    #         )
    #         vc[i][j] = 0.25 * (
    #             vf[2 * i][2 * j]
    #             + vf[2 * i + 1][2 * j]
    #             + vf[2 * i][2 * j + 1]
    #             + vf[2 * i + 1][2 * j + 1]
    #         )
    return uc, vc


def nonL(c_new, mu_new, nxt, nyt, dt, epsilon2, xright, xleft, yright, yleft, boundary):
    """
    NSO operator
    :param c_new: c at a time step
    :param mu_new: mu at a time step
    :param nxt: temp (number of grid points in x-direction, locally defined)
    :param nyt: temp (number of grid points in y-direction, locally defined)
    :return: ru, rw
    """
    ru = aux.dmatrix(nxt, nyt)
    rw = aux.dmatrix(nxt, nyt)

    lap_c = ll.laplace(c_new, nxt, nyt, xright, xleft, yright, yleft, boundary)
    lap_mu = ll.laplace(mu_new, nxt, nyt, xright,
                        xleft, yright, yleft, boundary)
    # for i in range(nxt):
    #     for j in range(nyt):
    #         ru[i][j] = c_new[i][j] / dt - lap_mu[i][j]
    #         rw[i][j] = mu_new[i][j] - \
    #             (c_new[i][j]) ** 3 + epsilon2 * lap_c[i][j]
    ru = c_new / dt - lap_mu
    rw = mu_new - (c_new) ** 3 + epsilon2 * lap_c
    return ru, rw


def defect(uf_new, wf_new, suf, swf, nxf, nyf, uc_new, wc_new, nxc, nyc, dt, epsilon2, xright, xleft, yright, yleft, boundary):
    ruc, rwc = nonL(uc_new, wc_new, nxc, nyc, dt, epsilon2,
                    xright, xleft, yright, yleft, boundary)
    ruf, rwf = nonL(uf_new, wf_new, nxf, nyf, dt, epsilon2,
                    xright, xleft, yright, yleft, boundary)
    ruf = suf - ruf
    rwf = swf - rwf
    rruf, rrwf = restrict_ch(ruf, rwf, nxc, nyc)
    duc = ruc + rruf
    dwc = rwc + rrwf
    return duc, dwc


def prolong_ch(uc, vc, nxc, nyc):
    uf = np.zeros((2 * nxc, 2 * nyc))
    vf = np.zeros((2 * nxc, 2 * nyc))
    for i in range(nxc):
        for j in range(nyc):
            uf[2 * i][2 * j] = uf[2 * i + 1][2 * j] = uf[2 * i][2 * j + 1] = uf[
                2 * i + 1
            ][2 * j + 1] = uc[i][j]
            vf[2 * i][2 * j] = vf[2 * i + 1][2 * j] = vf[2 * i][2 * j + 1] = vf[
                2 * i + 1
            ][2 * j + 1] = vc[i][j]
    return uf, vf


def vcycle(
    uf_new, wf_new, su, sw, nxf, nyf, ilevel, c_relax, xright, xleft, yright, yleft, dt, epsilon2, n_level, boundary,
):
    """
    FAS multigrid cycle
    """

    uf_new, wf_new = relax.relax(uf_new, wf_new, su, sw, nxf, nyf, c_relax,
                                 xright, xleft, yright, yleft, dt, epsilon2, boundary)
    if ilevel < n_level:
        nxc = int(nxf / 2)
        nyc = int(nyf / 2)
        uc_new, wc_new = restrict_ch(uf=uf_new, vf=wf_new, nxc=nxc, nyc=nyc)

        duc, dwc = defect(uf_new, wf_new, su, sw, nxf, nyf, uc_new, wc_new,
                          nxc, nyc, dt, epsilon2, xright, xleft, yright, yleft, boundary)

        uc_def = uc_new.copy()
        wc_def = wc_new.copy()

        uc_def, wc_def = vcycle(uc_def, wc_def, duc, dwc, nxc, nyc, ilevel + 1,
                                c_relax, xright, xleft, yright, yleft, dt, epsilon2, n_level, boundary)

        # uc_def = uc_def - uc_new
        uc_def -= uc_new
        wc_def -= wc_new

        uf_def, wf_def = prolong_ch(uc=uc_def, vc=wc_def, nxc=nxc, nyc=nyc)

        uf_new += uf_def
        wf_new += wf_def

        uf_new, wf_new = relax.relax(uf_new, wf_new, su, sw, nxf, nyf, c_relax,
                                     xright, xleft, yright, yleft, dt, epsilon2, boundary)

    return uf_new, wf_new


def cahn(c_old,
         c_new,
         mu,
         nx,
         ny,
         dt,
         solver_iter,
         tol,
         c_relax,
         xright,
         xleft,
         yright,
         yleft,
         epsilon2,
         n_level,
         boundary,
         suffix="",
         printres=True,
         pathname="",):
    it_mg2 = 0
    resid2 = 1
    sc, smu = source(c_old, nx, ny, dt, xright, xleft, yright, yleft, boundary)

    while it_mg2 < solver_iter and resid2 > tol:

        c_new, mu = vcycle(c_new, mu, sc, smu, nx, ny, 1, c_relax, xright,
                           xleft, yright, yleft, dt, epsilon2, n_level,
                           boundary)
        resid2 = error2.error2(c_old, c_new, mu, nx, ny, dt, xright, xleft,
                               yright, yleft, boundary)
        if printres:
            with open(f"{pathname}residual.csv", "a") as res:
                res.write(f"{resid2},\n")
        it_mg2 += 1

    return c_new
