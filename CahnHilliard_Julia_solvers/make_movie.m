outdir = "../output/output_julia";
n_relax = 4;
GridSize = 128;
dt = 1.00e-5;
max_it = 2000;

pathname = sprintf("%s/NMG_Julia_%d_dt_%.2e_Nx_%d_n_relax_%d_",outdir,max_it,dt, GridSize, n_relax);
t_out = readmatrix(strcat(pathname,"t_out.csv"));
t_out = transpose(t_out);
filename = strcat(pathname, "movie");
cd("/Users/smgroves/Documents/GitHub/CHsolvers_package/CahnHilliard_MATLAB_solvers")
ch_movie_from_file(strcat(pathname,"phi.csv"), t_out, GridSize,filename = filename);
