import numpy as np
cimport numpy as np
from . import aux_functions_NMG as aux
from . import laplacian as ll
from . import error2
from . import relax


cpdef source(c_old, int nx, int ny, double dt,
             double xright, double xleft, double yright, double yleft,
             str boundary):
    src_mu = np.zeros((nx, ny))
    ct = ll.laplace(c_old, nx, ny, xright, xleft, yright, yleft, boundary)
    src_c = np.asarray(c_old) / dt - np.asarray(ct)
    return src_c, src_mu


cpdef restrict_ch(uf, vf, int nxc, int nyc):
    uc = 0.25 * (np.asarray(uf)[::2, ::2] + np.asarray(uf)[1::2, ::2] + np.asarray(uf)[::2, 1::2] + np.asarray(uf)[1::2, 1::2])
    vc = 0.25 * (np.asarray(vf)[::2, ::2] + np.asarray(vf)[1::2, ::2] + np.asarray(vf)[::2, 1::2] + np.asarray(vf)[1::2, 1::2])
    return uc, vc


cpdef nonL(c_new, mu_new, int nxt, int nyt, double dt, double epsilon2,
           double xright, double xleft, double yright, double yleft,
           str boundary):
    lap_c = ll.laplace(c_new, nxt, nyt, xright, xleft, yright, yleft, boundary)
    lap_mu = ll.laplace(mu_new, nxt, nyt, xright, xleft, yright, yleft, boundary)
    ru = np.asarray(c_new) / dt - np.asarray(lap_mu)
    rw = np.asarray(mu_new) - (np.asarray(c_new)) ** 3 + epsilon2 * np.asarray(lap_c)
    return ru, rw


cpdef defect(uf_new, wf_new, suf, swf, int nxf, int nyf,
             uc_new, wc_new, int nxc, int nyc,
             double dt, double epsilon2,
             double xright, double xleft, double yright, double yleft,
             str boundary):
    ruc, rwc = nonL(uc_new, wc_new, nxc, nyc, dt, epsilon2,
                    xright, xleft, yright, yleft, boundary)
    ruf, rwf = nonL(uf_new, wf_new, nxf, nyf, dt, epsilon2,
                    xright, xleft, yright, yleft, boundary)
    ruf = np.asarray(suf) - np.asarray(ruf)
    rwf = np.asarray(swf) - np.asarray(rwf)

    rruf, rrwf = restrict_ch(ruf, rwf, nxc, nyc)
    duc = np.asarray(ruc) + np.asarray(rruf)
    dwc = np.asarray(rwc) + np.asarray(rrwf)
    return duc, dwc


cpdef prolong_ch(double[:, :] uc, double[:, :] vc, int nxc, int nyc):
    cdef int i, j
    cdef double[:, :] uf_view, vf_view
    uf = np.zeros((2 * nxc, 2 * nyc))
    vf = np.zeros((2 * nxc, 2 * nyc))
    uf_view = uf
    vf_view = vf
    for i in range(nxc):
        for j in range(nyc):
            uf_view[2*i, 2*j] = uf_view[2*i+1, 2*j] = uf_view[2*i, 2*j+1] = uf_view[2*i+1, 2*j+1] = uc[i, j]
            vf_view[2*i, 2*j] = vf_view[2*i+1, 2*j] = vf_view[2*i, 2*j+1] = vf_view[2*i+1, 2*j+1] = vc[i, j]
    return uf, vf


cpdef vcycle(double[:, :] uf_new, double[:, :] wf_new,
             double[:, :] su, double[:, :] sw,
             int nxf, int nyf, int ilevel, int c_relax,
             double xright, double xleft, double yright, double yleft,
             double dt, double epsilon2, int n_level, str boundary):
    cdef int nxc, nyc

    uf_new, wf_new = relax.relax(uf_new, wf_new, su, sw, nxf, nyf, c_relax,
                                  xright, xleft, yright, yleft, dt, epsilon2, boundary)
    if ilevel < n_level:
        nxc = nxf // 2
        nyc = nyf // 2
        uc_new, wc_new = restrict_ch(uf_new, wf_new, nxc, nyc)
        duc, dwc = defect(uf_new, wf_new, su, sw, nxf, nyf, uc_new, wc_new,
                          nxc, nyc, dt, epsilon2, xright, xleft, yright, yleft, boundary)
        uc_def = uc_new.copy()
        wc_def = wc_new.copy()
        uc_def, wc_def = vcycle(uc_def, wc_def, duc, dwc, nxc, nyc, ilevel + 1,
                                c_relax, xright, xleft, yright, yleft, dt, epsilon2, n_level, boundary)
        uc_def -= uc_new
        wc_def -= wc_new
        uf_def, wf_def = prolong_ch(uc_def, wc_def, nxc, nyc)
        uf_new += uf_def
        wf_new += wf_def
        uf_new, wf_new = relax.relax(uf_new, wf_new, su, sw, nxf, nyf, c_relax,
                                      xright, xleft, yright, yleft, dt, epsilon2, boundary)
    return uf_new, wf_new


cpdef cahn(double[:, :] c_old, double[:, :] c_new, double[:, :] mu,
           int nx, int ny, double dt, int solver_iter, double tol,
           int c_relax, double xright, double xleft, double yright, double yleft,
           double epsilon2, int n_level, str boundary,
           str suffix="", bint printres=True, str pathname=""):
    cdef int it_mg2 = 0
    cdef double resid2 = 1.0

    sc, smu = source(c_old, nx, ny, dt, xright, xleft, yright, yleft, boundary)
    while it_mg2 < solver_iter and resid2 > tol:
        c_new, mu = vcycle(c_new, mu, sc, smu, nx, ny, 1, c_relax, xright,
                           xleft, yright, yleft, dt, epsilon2, n_level, boundary)
        resid2 = error2.error2(c_old, c_new, mu, nx, ny, dt, xright, xleft,
                               yright, yleft, boundary)
        if printres:
            with open(f"{pathname}residual.csv", "a") as res:
                res.write(f"{resid2},\n")
        it_mg2 += 1
    return c_new