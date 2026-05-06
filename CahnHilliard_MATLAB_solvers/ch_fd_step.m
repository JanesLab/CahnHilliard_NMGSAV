function [u_new, info] = ch_fd_step(u_old, L, alpha, gamma, method, varargin)
%CH_FD_STEP One explicit or fully implicit Euler finite-difference step.
%
%   [u_new, info] = ch_fd_step(u_old, L, alpha, gamma, method, ...)
%
%   u_old is a column vector phi(:), L is the unscaled sparse FD Laplacian,
%   and alpha = dt*D. The implicit branch solves
%
%       phi^{n+1} - phi^n = alpha*L*mu^{n+1}
%       mu^{n+1} = (phi^{n+1})^3 - phi^{n+1} - gamma*L*phi^{n+1}
%
%   by damped Newton iteration.

p = inputParser;
addParameter(p, 'newton_tol', 1e-9);
addParameter(p, 'newton_max_iter', 20);
addParameter(p, 'linear_solver', 'backslash');
addParameter(p, 'line_search', true);
parse(p, varargin{:});
opt = p.Results;

method = lower(string(method));
info.iter = 0;
info.res = NaN;
info.ok = true;
info.method = char(method);

if method == "explicit" || method == "forward_euler"
    mu = u_old.^3 - u_old - gamma * (L * u_old);
    u_new = u_old + alpha * (L * mu);
    explicit_residual = u_new - u_old - alpha * (L * mu);
    info.res = norm(explicit_residual) / sqrt(numel(u_old));
    return;
end

if method ~= "implicit_euler" && method ~= "backward_euler" && method ~= "implicit"
    error('Unknown method: %s. Use explicit or implicit_euler.', method);
end

[u_new, info] = implicit_euler_step(u_old, L, alpha, gamma, ...
    opt.newton_tol, opt.newton_max_iter, lower(string(opt.linear_solver)), opt.line_search);
info.method = char(method);
end

function [u, info] = implicit_euler_step(u_old, L, alpha, gamma, tol, max_iter, linear_solver, do_line_search)
n = numel(u_old);
I = speye(n);
u = u_old;
info.iter = 0;
info.res = NaN;
info.ok = false;

for k = 1:max_iter
    [R, res] = residual(u, u_old, L, alpha, gamma);
    info.iter = k - 1;
    info.res = res;
    if res < tol
        info.ok = true;
        return;
    end

    dmu_du = spdiags(3*u.^2 - 1, 0, n, n) - gamma * L;
    J = I - alpha * (L * dmu_du);

    rhs = -R;
    if linear_solver == "gmres"
        [du, flag, relres] = gmres(J, rhs, [], 1e-10, min(200, n)); %#ok<ASGLU>
        if flag ~= 0
            warning('GMRES did not fully converge at Newton iter %d; flag=%d, relres=%g', k, flag, relres);
        end
    elseif linear_solver == "backslash"
        du = J \ rhs;
    else
        error('Unknown linear solver: %s. Use backslash or gmres.', linear_solver);
    end

    if do_line_search
        lambda = 1.0;
        accepted = false;
        while lambda > 1e-6
            u_trial = u + lambda * du;
            [~, res_trial] = residual(u_trial, u_old, L, alpha, gamma);
            if res_trial <= (1 - 1e-4 * lambda) * res || res_trial < tol
                u = u_trial;
                accepted = true;
                break;
            end
            lambda = 0.5 * lambda;
        end
        if ~accepted
            u = u + du;
        end
    else
        u = u + du;
    end
end

[~, res] = residual(u, u_old, L, alpha, gamma);
info.iter = max_iter;
info.res = res;
info.ok = (res < tol);
end

function [R, res] = residual(u, u_old, L, alpha, gamma)
mu = u.^3 - u - gamma * (L * u);
R = u - u_old - alpha * (L * mu);
res = norm(R) / sqrt(numel(u));
end
