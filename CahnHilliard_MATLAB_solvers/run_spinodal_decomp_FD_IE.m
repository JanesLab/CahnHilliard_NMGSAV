% FIGURE 1
indir = "./IC/";
outdir = "./output/output_MATLAB-periodic_FD_IE";
% make directories if they don't exist
if ~exist(outdir, 'dir')
    mkdir(outdir);
end
n_relax = 4;
m = 8;
GridSize = 128;
D = GridSize^2;
h = 1/GridSize;
epsilon = m * h/ (2 * sqrt(2) * atanh(0.9));
gamma = (m / (2 * sqrt(2) * atanh(0.9)))^2;
dt = 5.5e-6;
max_it = 2000;
boundary = 'periodic';
init_file = sprintf("%s/initial_phi_%d_smooth_n_relax_%d.csv",indir,GridSize, n_relax);
phi0 = readmatrix(init_file);
print_phi = true;
dt_out = 10;
ny = GridSize;
% #################################################
% RUN FD IE SOLVER 
% #################################################

pathname = sprintf("%s/FD_IE_MATLAB_%d_dt_%.2e_Nx_%d_n_relax_%d_dtout_1_",outdir,max_it,dt, GridSize, n_relax);
fprintf("Running FD IE solver with parameters: %s\n", pathname);
tStart_SAV = tic;
[phi_imp, st_imp] = ch_fd_solve(phi0, ...
    'method', 'implicit_euler', 'dt', dt, 't_iter', max_it, ...
    'D', D, 'gamma', gamma, 'boundary', boundary, 'save_every', dt_out, ...
    'newton_tol', 1e-9, 'newton_max_iter', 20, 'linear_solver', 'backslash');
elapsedTime = toc(tStart_SAV);

fid = fopen(sprintf('%s/Job_specs.csv', outdir), 'a+');
v = [string(datetime) "FD_IE_spinodal_decomp_smoothed_periodic_dtout_10_relaxation" "MATLAB" "FD_IE" GridSize epsilon dt max_it elapsedTime];
fprintf(fid, '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n', v);
fclose(fid);

% writematrix(phi_t(:,:,end),sprintf('%sfinal_phi.csv', pathname));
writematrix(st_imp.mass,sprintf('%smass.csv', pathname));
writematrix(st_imp.energy_norm,sprintf('%senergy.csv', pathname));

fprintf("Creating movie\n");
filename = strcat(pathname, "movie");

ch_movie(phi_imp,st_imp.t, filename = filename);



