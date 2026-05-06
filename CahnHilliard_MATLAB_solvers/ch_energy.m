function E = ch_energy(phi, gamma, boundary)
%CH_ENERGY Discrete Cahn-Hilliard energy used for FD diagnostics.
%
%   E = sum 1/4*(phi^2 - 1)^2 + gamma/2 * sum |grad_h phi|^2
%
%   The gradient term uses the same unscaled mesh convention as the current
%   finite-difference code. For Neumann boundaries, only interior forward
%   differences are counted. For periodic boundaries, wrapped differences are
%   included.

if nargin < 3 || isempty(boundary)
    boundary = 'neumann';
end
boundary = lower(string(boundary));

f = 0.25 * (phi.^2 - 1).^2;
E_local = sum(f(:));

if boundary == "periodic"
    dx = circshift(phi, -1, 1) - phi;
    dy = circshift(phi, -1, 2) - phi;
    E_grad = sum(dx(:).^2) + sum(dy(:).^2);
elseif boundary == "neumann"
    dx = diff(phi, 1, 1);
    dy = diff(phi, 1, 2);
    E_grad = sum(dx(:).^2) + sum(dy(:).^2);
else
    error('Unknown boundary condition: %s. Use neumann or periodic.', boundary);
end

E = E_local + 0.5 * gamma * E_grad;
end
