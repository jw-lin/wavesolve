from wavesolve.FEsolver import solve_waveguide_vec
# default Nmax is 6
eigvals_ord1,eigvecs_ord1 = solve_waveguide_vec(mesh_ord1,wl,IOR_dict,Nmax=7)