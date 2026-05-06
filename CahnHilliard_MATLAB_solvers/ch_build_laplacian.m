function L = ch_build_laplacian(nx, ny, boundary)
%CH_BUILD_LAPLACIAN Sparse 2-D unscaled finite-difference Laplacian.
%
%   L = ch_build_laplacian(nx, ny, boundary)
%
%   The returned matrix follows the same mesh convention as the original
%   finite-difference code: the mesh spacing is not included in L. Scaling by
%   h^{-2} is supplied separately through D = GridSize^2 when the domain
%   length is one.
%
%   boundary = 'neumann'  -> zero normal derivative, matching laplacian.m
%   boundary = 'periodic' -> wrap-around boundary
%
%   Vectorization follows MATLAB column-major order: phi(:).

if nargin < 3 || isempty(boundary)
    boundary = 'neumann';
end
boundary = lower(string(boundary));

Tx = one_dim_laplacian(nx, boundary);
Ty = one_dim_laplacian(ny, boundary);

L = kron(speye(ny), Tx) + kron(Ty, speye(nx));
end

function T = one_dim_laplacian(n, boundary)
e = ones(n, 1);
T = spdiags([e -2*e e], -1:1, n, n);

if boundary == "periodic"
    T(1, n) = 1;
    T(n, 1) = 1;
elseif boundary == "neumann"
    T(1, 1) = -1;
    T(n, n) = -1;
else
    error('Unknown boundary condition: %s. Use neumann or periodic.', boundary);
end
end
