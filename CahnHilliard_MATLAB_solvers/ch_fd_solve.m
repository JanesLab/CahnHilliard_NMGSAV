function [phi_t, stats] = ch_fd_solve(phi0, varargin)
%CH_FD_SOLVE Explicit or fully implicit Euler FD solver for Cahn-Hilliard.
%
%   [phi_t, stats] = ch_fd_solve(phi0, 'method', 'implicit_euler', ...)
%
%   Supported time discretizations:
%       method='explicit'        : forward Euler, matching the legacy FD code
%       method='implicit_euler'  : fully implicit backward Euler
%
%   stats.energy is unnormalized; stats.energy_norm = stats.energy/E(0).

[nx, ny] = size(phi0);
default_gamma = (8 / (2 * sqrt(2) * atanh(0.9)))^2;

p = inputParser;
addParameter(p, 'method', 'explicit');
addParameter(p, 'dt', 5.5e-6);
addParameter(p, 't_iter', 1000);
addParameter(p, 'D', nx^2);
addParameter(p, 'gamma', default_gamma);
addParameter(p, 'boundary', 'neumann');
addParameter(p, 'save_every', 10);
addParameter(p, 'newton_tol', 1e-9);
addParameter(p, 'newton_max_iter', 20);
addParameter(p, 'linear_solver', 'backslash');
addParameter(p, 'line_search', true);
addParameter(p, 'stop_threshold', 1e8);
addParameter(p, 'verbose', false);
parse(p, varargin{:});
opt = p.Results;

method = lower(string(opt.method));
boundary = lower(string(opt.boundary));
save_every = max(1, round(opt.save_every));
L = ch_build_laplacian(nx, ny, boundary);
alpha = opt.dt * opt.D;

u = phi0(:);

n_save = floor(opt.t_iter / save_every) + 1;
if mod(opt.t_iter, save_every) ~= 0
    n_save = n_save + 1;
end
phi_t = zeros(nx, ny, n_save);
phi_t(:, :, 1) = phi0;

stats = initialize_stats(opt.t_iter, opt.dt, u, phi0, opt.gamma, boundary, save_every);

save_idx = 1;
tic;
for step = 1:opt.t_iter
    [u, step_info] = ch_fd_step(u, L, alpha, opt.gamma, method, ...
        'newton_tol', opt.newton_tol, ...
        'newton_max_iter', opt.newton_max_iter, ...
        'linear_solver', opt.linear_solver, ...
        'line_search', opt.line_search);

    phi = reshape(u, nx, ny);
    stats.t(step + 1) = step * opt.dt;
    stats.mass(step + 1) = mean(u);
    stats.energy(step + 1) = ch_energy(phi, opt.gamma, boundary);
    stats.newton_iter(step) = step_info.iter;
    stats.newton_res(step) = step_info.res;
    stats.newton_ok(step) = step_info.ok;

    if ~all(isfinite(u)) || max(abs(u)) > opt.stop_threshold
        stats.blowup_step = step;
        if opt.verbose
            fprintf('Stopped at step %d due to nonfinite or huge phi.\n', step);
        end
        stats = trim_stats(stats, step);
        break;
    end

    if mod(step, save_every) == 0 || step == opt.t_iter
        save_idx = save_idx + 1;
        phi_t(:, :, save_idx) = phi;
        stats.saved_steps(save_idx) = step;
        stats.saved_t(save_idx) = step * opt.dt;
    end

    if opt.verbose && mod(step, max(1, round(opt.t_iter/10))) == 0
        fprintf('%s step %d/%d, E/E0 = %.6g, mass err = %.3e\n', ...
            method, step, opt.t_iter, stats.energy(step+1)/stats.energy(1), ...
            stats.mass(step+1)-stats.mass(1));
    end
end

stats.elapsed = toc;
phi_t = phi_t(:, :, 1:save_idx);
stats.saved_steps = stats.saved_steps(1:save_idx);
stats.saved_t = stats.saved_t(1:save_idx);
stats.energy_norm = stats.energy / stats.energy(1);
end

function stats = initialize_stats(t_iter, dt, u, phi0, gamma, boundary, save_every)
stats.t = zeros(t_iter + 1, 1);
stats.mass = zeros(t_iter + 1, 1);
stats.energy = zeros(t_iter + 1, 1);
stats.energy_norm = zeros(t_iter + 1, 1);
stats.newton_iter = nan(t_iter, 1);
stats.newton_res = nan(t_iter, 1);
stats.newton_ok = true(t_iter, 1);
stats.blowup_step = NaN;
stats.elapsed = NaN;
stats.save_every = save_every;
stats.saved_steps = nan(floor(t_iter/save_every) + 2, 1);
stats.saved_t = nan(floor(t_iter/save_every) + 2, 1);
stats.saved_steps(1) = 0;
stats.saved_t(1) = 0;
stats.t(1) = 0;
stats.mass(1) = mean(u);
stats.energy(1) = ch_energy(phi0, gamma, boundary);
end

function stats = trim_stats(stats, step)
stats.t = stats.t(1:step+1);
stats.mass = stats.mass(1:step+1);
stats.energy = stats.energy(1:step+1);
stats.newton_iter = stats.newton_iter(1:step);
stats.newton_res = stats.newton_res(1:step);
stats.newton_ok = stats.newton_ok(1:step);
end
