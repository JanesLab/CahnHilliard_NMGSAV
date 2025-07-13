function relax!(c_new, mu_new, su, sw, nxt, nyt, c_relax, xright, xleft, yright, yleft, dt, epsilon2, boundary)
    ht2 = ((xright - xleft) / nxt)^2
    a = MVector{4,Float64}(undef)
    f = MVector{2,Float64}(undef)
    for iter in 1:c_relax
        for i in 1:nxt
            for j in 1:nyt
                if boundary == "neumann"
                    if i > 1 && i < nxt
                        x_fac = 2.0
                    else
                        x_fac = 1.0
                    end
                    if j > 1 && j < nyt
                        y_fac = 2.0
                    else
                        y_fac = 1.0
                    end
                elseif boundary == "periodic"
                    x_fac = 2.0
                    y_fac = 2.0
                end
                a[1] = 1 / dt
                a[2] = (x_fac + y_fac) / ht2
                a[3] = -(x_fac + y_fac) * epsilon2 / ht2 - 3 * (c_new[i, j])^2
                a[4] = 1.0

                f[1] = su[i, j]
                f[2] = sw[i, j] - 2 * (c_new[i, j])^3

                if i > 1
                    f[1] += mu_new[i-1, j] / ht2
                    f[2] -= epsilon2 * c_new[i-1, j] / ht2
                elseif boundary == "periodic"
                    f[1] += mu_new[nxt, j] / ht2 #changed
                    f[2] -= epsilon2 * c_new[nxt, j] / ht2 #changed
                end
                if i < nxt
                    f[1] += mu_new[i+1, j] / ht2
                    f[2] -= epsilon2 * c_new[i+1, j] / ht2
                elseif boundary == "periodic"
                    f[1] += mu_new[1, j] / ht2 #changed
                    f[2] -= epsilon2 * c_new[1, j] / ht2 #changed
                end
                if j > 1
                    f[1] += mu_new[i, j-1] / ht2
                    f[2] -= epsilon2 * c_new[i, j-1] / ht2
                elseif boundary == "periodic"
                    f[1] += mu_new[i, nyt] / ht2 #changed
                    f[2] -= epsilon2 * c_new[i, nyt] / ht2 #changed
                end
                if j < nyt
                    f[1] += mu_new[i, j+1] / ht2
                    f[2] -= epsilon2 * c_new[i, j+1] / ht2
                elseif boundary == "periodic"
                    f[1] += mu_new[i, 1] / ht2 #changed
                    f[2] -= epsilon2 * c_new[i, 1] / ht2 #changed
                end
                det = a[1] * a[4] - a[2] * a[3]
                c_new[i, j] = (a[4] * f[1] - a[2] * f[2]) / det
                mu_new[i, j] = (-a[3] * f[1] + a[1] * f[2]) / det

            end
        end
    end
    # return c_new, mu_new
end
