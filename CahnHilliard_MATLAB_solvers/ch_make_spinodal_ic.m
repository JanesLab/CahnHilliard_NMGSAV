function phi0 = ch_make_spinodal_ic(GridSize, plus_fraction, varargin)
%CH_MAKE_SPINODAL_IC Deterministic spinodal initial condition for tests.
%
%   phi0 = ch_make_spinodal_ic(GridSize, plus_fraction, ...)
%
%   Name-value options:
%       'seed'         : RNG seed, default 1234
%       'boundary'     : 'neumann' or 'periodic', default 'neumann'
%       'n_relax'      : number of heat-smoothing passes, default 4
%       'smooth_theta' : heat-smoothing strength, default 0.1
%
%   This is deliberately self-contained so the reviewer tests do not depend
%   on absolute paths to initial-condition CSV files.

p = inputParser;
addParameter(p, 'seed', 1234);
addParameter(p, 'boundary', 'neumann');
addParameter(p, 'n_relax', 4);
addParameter(p, 'smooth_theta', 0.1);
parse(p, varargin{:});
opt = p.Results;

rng(opt.seed);
N = GridSize * GridSize;
phi0 = -ones(GridSize, GridSize);
idx = randperm(N, round(plus_fraction * N));
phi0(idx) = 1;

if opt.n_relax > 0
    L = ch_build_laplacian(GridSize, GridSize, opt.boundary);
    u = phi0(:);
    for k = 1:opt.n_relax
        u = u + opt.smooth_theta * (L * u);
    end
    phi0 = reshape(u, GridSize, GridSize);
end
end
