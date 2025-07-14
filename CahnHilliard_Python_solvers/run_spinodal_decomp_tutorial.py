import CHsolvers as ch
import pandas as pd
import os
from datetime import datetime
import numpy as np


indir = "./IC/"
boundary = "periodic"

outdir = f"./output/output_python-{boundary}/"

method = "SAV"
n_relax = 4
GridSize = 128
h = 1/GridSize
m = 8
epsilon = m * h / (2 * np.sqrt(2) * np.arctanh(0.9))

dt = 5.5E-06
max_it = 20
printphi = True
dt_out = max_it  # output every 10 timesteps


# phi0 = ch.init.initialization_from_file(f"{indir}initial_phi_$(GridSize)_smooth_n_relax_{n_relax}.csv",
#                                         GridSize, GridSize, delim=",", transpose_matrix=False)

phi0 = np.loadtxt(
    f"{indir}initial_phi_{GridSize}_smooth_n_relax_{n_relax}.csv", delimiter=",")

pathname = f"{outdir}{method}_Python_{max_it}_dt_{dt:.2e}_Nx_{GridSize}_n_relax_{n_relax}_{boundary}_"

date_time = datetime.now().strftime("%m/%d/%Y %H:%M:%S")

# note that if using the time_and_mem decorator (in CahnHilliard_NMG), it will return a results dictionary.
# TODO remove decorator for final package
results_dict = ch.SAV.CahnHilliard_SAV(phi0, t_iter=max_it, dt=dt, dt_out=dt_out, m=m,
                                       boundary=boundary, printphi=printphi, printres=True, pathname=pathname)

# t_out, phi_t, delta_mass_t, E_t
# rename results_dict keys from results1, results2, etc. to t_out, phi_t, delta_mass_t, E_t
results_dict["t_out"] = results_dict.pop("results1")
results_dict["phi_t"] = results_dict.pop("results2")
results_dict["delta_mass_t"] = results_dict.pop("results3")
results_dict["E_t"] = results_dict.pop("results4")

# Save results
np.savetxt(f"{pathname}t_out.csv", results_dict["t_out"], delimiter=",")
np.savetxt(f"{pathname}delta_mass_t.csv",
           results_dict["delta_mass_t"], delimiter=",")
np.savetxt(f"{pathname}E_t.csv", results_dict["E_t"], delimiter=",")

# Save phi_t by appending to file
# if printphi == False:
#     with open(f"{pathname}phi.csv", "w") as f:
#         for t in range(max_it+1):
#             for i in range(GridSize):
#                 for j in range(GridSize):
#                     f.write(f"{results_dict['phi_t'][i,j,t]},")
#                 f.write("\n")  # end of each row


T = {}
T['date'] = date_time
T['name'] = "spinodal_decomp_smoothed_save_variable_dtout_20_print"
T['language'] = "Python"
T['solver'] = method
T['nx'] = GridSize
T['epsilon'] = epsilon
T['dt'] = dt
T['tol'] = np.nan
T['timesteps'] = max_it
T['max_it_CH'] = np.nan
T['output_name'] = pathname
T['time (secs)'] = results_dict["computation_time_Sec"]
T['mem_allocated (MB)'] = results_dict["memory_usage_MB"]
T["bc"] = boundary

T = pd.DataFrame([T])
file = "./Job_specs.csv"
if not os.path.isfile(file):
    T.to_csv(file)
else:
    with open(file, "a") as f:
        T.to_csv(f, header=False, index=False)
