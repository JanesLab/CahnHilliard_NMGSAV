
# laplacian function: laplacian(m, nx, ny, h2)
function laplace!(lap_a, a, nxt, nyt, xright, xleft, yright, yleft, boundary)

    ht2 = ((xright - xleft) / nxt)^2
    for i in 1:nxt
        for j in 1:nyt
            if i > 1
                dadx_L = (a[i, j] - a[i-1, j])
            else
                if boundary == "neumann"
                    dadx_L = 0
                elseif boundary == "periodic"
                    dadx_L = a[i, j] - a[nxt, j] #changed
                end
            end
            if i < nxt
                dadx_R = (a[i+1, j] - a[i, j])
            else
                if boundary == "neumann"
                    dadx_R = 0
                elseif boundary == "periodic"
                    dadx_R = a[1, j] - a[i, j] #changed
                end
            end
            if j > 1
                dady_B = (a[i, j] - a[i, j-1])
            else
                if boundary == "neumann"
                    dady_B = 0
                elseif boundary == "periodic"
                    dady_B = a[i, j] - a[i, nyt] #changed
                end
            end
            if j < nyt
                dady_T = (a[i, j+1] - a[i, j])
            else
                if boundary == "neumann"
                    dady_T = 0
                elseif boundary == "periodic"
                    dady_T = a[i, 1] - a[i, j] #changed
                end
            end
            lap_a[i, j] = (dadx_R - dadx_L + dady_T - dady_B) / ht2
        end
    end
    return lap_a
end
